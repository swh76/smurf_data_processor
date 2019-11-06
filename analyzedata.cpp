// extracts and analyzes data from a smurf output file to a text file.
// can operate with very large files

// input switches
// -I  input file name
// -O output file name
// -L channel list file name. This is a text file with a list of channels to read
// -A number of averages / downsample, simple rectangular average, not anti-aliased
// -D set to 1 to turn on diagnostic timing,

// output columnns
// 1st column  Time in seconds, small, offset from 0. this is from the internal SMURF timing system
// If diagnostics is off, 2nd through last columns are data, floating point, space separated
// If diagnoistics is on -D 1,   then:
// 2nd column,  Internal SMURF frame counter - eg flux ramp counter
// 3rd column,  MCE syncbox counter (if used, otherwise 0, or some other instrument time)
// 4th column, Unix time in nanoseconds when data was taken,  Divide by 1e9 to get time since 1970.  
// 5th and additional columns data in order given in list file. 

// for multi part files, can ccat all filess together      cat datafile.dat.00* > datafile.dat
// Then read the giant resulting file
// note, using too many channels at once will generate insanely huge text files


// example command
// ./analyzedata -I inputfilename.dat -O oututfilename.txt -L channellistfile.txt -A 10 -D 1
// that read inputfilename.dat,  produces outpufilename.txt   using the channelllistfile.txt  and turns on diagnostics.
// output is averaged and downsampled by 10

// build with  g++ -o analyzedata analyzedata.cpp
// build only takes a second



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
typedef uint32_t uint;


const uint buffsize = 100000000; // arbitrary for now.
const uint header_length = 128; // nubmer of bytes in header
const uint header_channel_offset = 4;  // bytes
const uint time_lower_offset = 72; // offset to 32 bit timing sysetm counter nanoseconds
const uint time_upper_offset = 76;  // offset to 32 bit timing system contuer seconds (maybe broken)
const uint mce_counter_offset = 96;
const uint frame_counter_offset = 84;
const uint data_size = 4; // data size bytes
const uint max_channels = 4096; // will really use 528
const uint unix_time_offset = 48; // 64 bit number

int exit_error( const char *errstring);

int exit_error(const char* errstring)
{
  printf("exiting, %s", errstring);
  exit(0);
}


int main(int argc, char* argv[])
{
	char infilename[256];
	char outfilename[256];
	char listfilename[256]; // name of mask file
	uint chlist[max_channels]; 
	char *tmpname; // simplifies some notiation
	//uint first_channel;
	//uint last_channel;
	//uint used_channels; 
	uint averages;
	uint read_size; // size of single read frame
	uint read_block; // number of bytes to read on each read operation.
	int fdin; // File descriptor input
	FILE* fpout; // File pointer output
	FILE* fplist; // file pointer for list file
	uint8_t* buffer;
	uint32_t num_channels;  // from header, number of channels in data
	uint output_channels;  // number of channels to write
	uint read_frames; // number of frames to read each call
	uint proc_frames;
	uint bytes_read;
	uint8_t* bufptr;
	int32_t data_out[max_channels]; // output data 
	uint32_t lower_time_counter;
	uint32_t upper_time_counter;
	uint32_t initial_upper_time;
	uint32_t initial_lower_time;
	uint ctr;
	uint tmp;
	uint avgcnt; // counter for average cycle
	uint outn; // output number
	double x, dtu, dtl;  // used for output, to simply syntax
	uint diagmode;
	uint32_t mcecounter;
	uint32_t framecounter;
	uint64_t unix_time;

	buffer = (uint8_t*)malloc(buffsize * sizeof(uint8_t));
	strcpy(infilename, "data.dat");
	strcpy(outfilename, "outfile.dat");
	strcpy(listfilename, "listfile.txt"); 

	averages = 1;
	diagmode = 0;  // normal mode 


	for (ctr = 1; ctr < argc; ctr++) // loop over input arguments
	  { 
	    tmpname = argv[ctr]; // just make code clearer
	    if((tmpname[0] == '-') && (argc > ctr))  // found an identfier and have another input
	      {
		switch(tmpname[1])   // what input control type
		  {
		  case 'I' :
		    strcpy(infilename, argv[++ctr]);
		    printf("found input file identifier, name is %s\n", infilename);
		    break;
		  case 'O' :
		    strcpy(outfilename, argv[++ctr]);
		    printf("found output file identifier, name is %s\n", outfilename);
		    break;
		  case 'L' :
		    strcpy(listfilename, argv[++ctr]);
		    printf("found channel list file identifier, name is %s\n", listfilename);
		    break;
		  case 'A':  // number of averages
		    averages = atoi(argv[++ctr]);
		    if(averages < 1) averages = 1; // no negative averages
		    printf("averages changed to changed to %d \n", diagmode); 
		    break;
		  case 'D':  // diagnostic mode
		    diagmode = atoi(argv[++ctr]);
		    printf("diagnostic mode changed to %d \n", diagmode); 
		    break;
		default:
		  printf("default found"); 
		  } 
	      }
	  }





	if (0 >= (fdin = open(infilename, O_RDONLY))){ printf("couldn not open input file %s \n", infilename);  exit(0);}
	if (!(fpout = fopen(outfilename, "w"))) {printf("couldn not open output file %s \n", outfilename);  exit(0); }	
	if (!(fplist = fopen(listfilename, "r"))) {printf("couldn not open list file %s \n", listfilename);  exit(0);}
	output_channels = 0;
	memset(chlist, 0, sizeof(uint) * max_channels); // clear channel list
        for(ctr = 0; ctr < max_channels; ctr++)
	  {
	    if (EOF == fscanf(fplist, "%u", chlist+ctr)) break;  // end of file or similar break, otherwise fill in value
	  }
	fclose(fplist);  // close file
	output_channels = ctr; // how many are we usign






	read(fdin, buffer, header_length); // read first header (will assume all ar the same
	close(fdin); // ugly way to rewind -must be a better way

	num_channels = *((uint32_t*)(buffer + header_channel_offset)); // sets read block size
	initial_upper_time = *((uint32_t*)(buffer + time_upper_offset)); // 64 bit counter 
	initial_lower_time = *((uint32_t*)(buffer + time_lower_offset));

	read_size = num_channels * data_size + header_length;  // all in bytes
	printf("arc = %d, fdin = %d, infile = %s, outfile = %s, numchans = %u\n", argc, fdin, infilename, outfilename, num_channels);

	tmp = buffsize / read_size / averages;
	read_frames = tmp * averages;  // number of frames at each read
	read_block = read_frames * read_size; // total bytes to read each time




	if (!(fdin = open(infilename, O_RDONLY))) return(0); // now open again. (really dumb)
	outn = 0;  // which output sample are we on.
	for (ctr = 0; ctr < 1000000000; ctr++)           //MAIN READ LOOP
	{
		bytes_read = read(fdin, buffer, read_block); // read-in block of data to process
		proc_frames = bytes_read / read_size;
		printf("read block = %u, bytes read = %u , proc_frames = %u\n", ctr, bytes_read, proc_frames);
		avgcnt = 0;
		for (uint j = 0; j < proc_frames; j++)
		{
			bufptr = buffer + j * read_size; // pointer to current buffer

			//for (uint n = first_channel; n <= last_channel; n++)
			for (uint m = 0; m < output_channels; m++)
			{
			  uint n = chlist[m];  // channel look up
			  if((n < 0) || (n > max_channels)) n = 0; // set limits
			  data_out[n] += *((int32_t*)(buffer + header_length + 4 * n + read_size * j)); // sum data	
			}
			avgcnt++;  // increment average count
			if (avgcnt == averages) // done averaging, time to write
			{
				unix_time = *((uint64_t*)(buffer + unix_time_offset + read_size * j));
				upper_time_counter = *((uint32_t*)(buffer + time_upper_offset + read_size * j)); // 32 bit counter
				lower_time_counter = *((uint32_t*)(buffer + time_lower_offset + read_size * j));
				dtl = (double)(lower_time_counter);
				dtu = (double)(upper_time_counter)-(double)(initial_upper_time);  // now a double
				framecounter = *((uint32_t*)(buffer + frame_counter_offset + read_size * j));
				mcecounter = *((uint32_t*)(buffer + mce_counter_offset + read_size * j));
				fprintf(fpout, "%12.6f  ", dtu + dtl / 1e9); // REAL
				if (diagmode == 1)
				{
					fprintf(fpout, "%10d   %10d   %lu ", framecounter, mcecounter, unix_time);
				}
				//for (uint n = first_channel; n <= last_channel; n++)
				for (uint m = 0; m < output_channels; m++)
				{
				  uint n = chlist[m];  // channel look up
				  if((n < 0) || (n > max_channels)) n = 0; // set limits
				  x = data_out[n]; // convert to float
				  fprintf(fpout, "%12.2f ", x / (double)averages);
				  data_out[n] = 0;  // clear
				}
				fprintf(fpout, "\n"); // new line
				avgcnt = 0;
				outn = 0;
			}
		}
		if (bytes_read < read_block) break; // done
	}

	if (fdin) close(fdin);
	if (fpout) fclose(fpout);


}
