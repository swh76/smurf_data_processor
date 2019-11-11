# File reader for smurf data analysis
# should be replaced with pysmurf reader
"""file_reader takes varous input smurf files.
Output is an array [data point, channel]
column 0:  timing
column 1-3:   if R['diagtiming'] is set, these are other timing diagnostics, if not these are data
additional columns: data

R['temp_directory'] usually set to /tmp   stores intermediate files
R['data_directory']  directory where data is stored
R['run_name']  files are assumed to be in the form of  abcde.dat,    or abcde.dat.part_0xxxx  
where abcde is R['run_name']


FILE EXTRACTION (separate C code)
if R'cat_files'] is True,  then a system call is used to cat all files of the form abcde.dat.part_0xxxx into 
a file abdce.dat  in the temporay directory

R['num_averages'] uses the data extraction program to dowsample and average data by the listed ratio

R['list_file'] is the name of the file listing which channels are toe be extracted

R['extract_program'] is the executable for the C file reader.  usuall ./analyzedata  
This is a C program that reads the (concatinated) smurf .dat file and writes text files out
This is an updated version of extractdata.cpp

R['textfilename'] is the text file for the extracted data

R['diagtiming'] = True causes 3 extra timing diagnostic outputs to be written. 

The help file for the extract program follows:

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
// If diagnoistics is on (use switch     -D 1 ),   then:
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
// throughput 1GB/s for input files on SSDs



"""
import subprocess   # used for the cat command
import numpy as np
import glob
import time
from SmurfFileReader import SmurfStreamReader

def new_file_reader(R):

    stime = time.time()

    R['txtfname'] = R['temp_directory'] +  R['run_name'] + '.txt' # extract to temp diretory

    if R['extract_txt']:
        chanList = []

        if R['num_averages'] == 0:
            numAverages = 1
        else:
            numAverages = R['num_averages']

        # Read in channel list
        with open(R['list_file']) as f:
            for line in f.readlines():
                if line.rstrip().isdigit():
                    chanList.append(int(line))

        # Init channel sums
        initialTime = None

        if R['cat_files']:
            fileList = sorted(glob.glob(R['data_directory'] + R['run_name'] + '.dat.part_0*'))
        else:
            fileList = [R['data_directory'] + R['run_name'] + '.dat']

        with open(R['txtfname'],'w') as of, SmurfStreamReader(fileList,isRogue=R['rogue_format']) as fr:

            for header,data in fr.records():

                if initialTime is None:
                    initialTime = header.counter_2
                    chanSum = np.zeros(len(data),np.float)
                    avgCount = 0

                if numAverages != 1:
                    chanSum += data

                avgCount += 1

                if avgCount == numAverages:

                    if numAverages != 1:
                        chanAvg = chanSum / avgCount
                        chanSum = np.zeros(len(data),np.float)
                    else:
                        chanAvg = data

                    avgCount = 0

                    # Time
                    diff = header.counter_2 - initialTime
                    real = (diff >> 32) + (diff & 0xFFFFFFFF) / 1e9
                    line = '{0:.6f}'.format(real)

                    if R['diag_timing']:
                        line += ' {} {} {}'.format(header.frame_counter,header.external_time,header.timestamp)

                    for c in chanList:
                        line += ' {0:.2f}'.format(chanAvg[c])

                    of.write(line + '\n')

    D = np.loadtxt(R['txtfname']) # numpy import of text file

    print("new_file_reader completed in {} seconds".format(time.time() - stime))
    return(D)

def file_reader(R):

    stime = time.time()

    if R['cat_files']:
        datfname = R['temp_directory'] +  R['run_name'] + '.dat'
        catcmd = 'cat ' + R['data_directory']+ R['run_name'] + '.dat.part_0*' + ' > ' + datfname
        print('catinating files using:  ' + catcmd) 
        x = subprocess.call(catcmd, shell=True)  # horrible kludge to combine files into one
        print('done cat files')
    else:
        datfname = R['data_directory'] +  R['run_name'] + '.dat'
    R['txtfname'] = R['temp_directory'] +  R['run_name'] + '.txt' # extract to temp diretory
    if R['extract_txt']:
        print(R['compile_command'])
        compile_comand = R['compile_command']
        print("compile = ")
        print(compile_comand)
        x = subprocess.call(compile_comand, shell=True)
        print("compile done")
        x = subprocess.call('chmod 777 ' + R['extract_program'], shell=True)   # super ugly but need to execute after compile
        extractcmd = R['extract_program'] + ' -I ' + datfname + ' -O ' + R['txtfname'] + ' -L ' + R['list_file']
        if R['diag_timing']:
            extractcmd = extractcmd + ' -D 1 '  # include diagnostics file, (that is the number 1)
        if R['num_averages'] > 1:
            extractcmd = extractcmd + ' -A ' + str(R['num_averages'])
            print('extracting: ' +  extractcmd)
        x = subprocess.call(extractcmd, shell=True) 
    D = np.loadtxt(R['txtfname']) # numpy import of text file

    print("file_reader completed in {} seconds".format(time.time() - stime))
    return(D)

class channel_mapping:
    def __init__(self, maskfilename, listfilename):
        self.map =np.loadtxt(maskfilename) # loads list
        self.list = np.loadtxt(listfilename)
        self.use_gcpmap = False
        
    def get_channel(self, n):
        listn = int(self.list[n])
        smurf_chan = int(self.map[listn]) # look up what smurf channel
        band = int(smurf_chan / 512)
        bchan = int(smurf_chan % 512)
        return( smurf_chan, band, bchan)
        
    def use_gcp_map(self, mapfilename): #map gcp 
        self.gcpmap = np.loadtxt(mapfilename)
        self.use_gcpmap = True 
        
    def gcp_index(self, n):
        print("gcp n = ", n)
        nz = int(self.list[n])
        print("gcp nz = ", nz)
        print("gcp map size = ", self.map.shape)
        mz = int(self.map[nz])
        print("gcp mz = ", mz)
        ns = np.where(self.gcpmap[:,1] == mz)
        nout = 0
        print("gcp ns = ", ns[0])
        print("length of ns = ", len(ns[0]))
        if len(ns[0] > 0):
            nx = int(ns[0]) # find index
            print("gcp nx = ", nx)
            nout = int(self.gcpmap[nx,0])
            print("gcp out = ", nout)
        return(nout)  # return value
        
            
