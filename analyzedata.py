#analyzedata.py
 
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import datetime
import file_reader

max_pairs = 100 # just to limit for now
 
    

        
    
def run(D, R): # R is a dictinoary that has control parameters

    run_string = R['run_name']
    show_plots = R['show_plots']
 
    filename = R['txtfname']
    maskfilename = R['mask_file'] 
    listfilename = R['list_file']
    reportfilename = R['report_file']
    first_sample = 0
    last_sample = 0  # zero lenth means use internal
 
    
    samples = last_sample - first_sample
    enable_spectrum_plot = True
    enable_time_plot = True
    welch_blocks = R['welch_blocks']
    timing_channels = 4  # number of channels used for timing
    reportfile = open(reportfilename, 'w' )  # for output info
    
    description = R['description']
    print('in analyze data descripton: ', description)
        
    
    
    reportfile.write("run name = " + run_string + '\n')
    reportfile.write(description + '\n')
    
    
    
    plt.close('all')
    Cmap = file_reader.channel_mapping(maskfilename, listfilename) # create channel mapping class
    Cmap.use_gcp_map('GCP_map.txt')

    ds = D.shape  #data shape
    if R['diag_timing']:
        nt = 4
    else:
        nt = 1
    if R['fake_data']:
        if R['fake_data'] >= 1: # mode 1 fake data, random         
            tottime = D[-1, 0] - D[0,0] # time difference
            print("**********************************************")
            print("total time = ", tottime)
            totsamples = ds[0]
            print("total sample = ", totsamples)
            D[:,nt:ds[1]] = R['fake_data'] * np.random.randn(ds[0], ds[1]-nt) * np.sqrt(.5*totsamples / tottime)    # overwrite iwth fake data      
        
    original_samples = ds[0]
    print("original samples = ", original_samples)
    chopfront = int(R['chop_front'] * original_samples)
    chopend = int(R['chop_back'] * original_samples)
    D = D[chopfront : original_samples-chopend,:]
    ds = D.shape  #data shape
    num_pairs = int((ds[1]-4)/2)
    total_samples = ds[0]    
    if samples == 0 :
        samples = total_samples
        first_sample = 0
        last_sample = total_samples
    R['total_samples'] = total_samples # add ths for scaling later 
    print("total samples = ", total_samples);
    A = np.zeros([samples,num_pairs])
    B = np.zeros([samples,num_pairs])
    Ma = np.zeros(num_pairs) # offsets
    Mb = np.zeros(num_pairs)
    print("pairs = ", num_pairs)
    T = D[first_sample:last_sample, 0] # time
    data_time = D[0,3]/ 1.0e9
    dtmp = datetime.datetime.utcfromtimestamp(int(data_time))
    time_string = dtmp.strftime('%c')
    time_string = " " + time_string + " " 
    
    
    # setup title
    title = run_string + " " + time_string + "\n" + description
    
    print("time = ", time_string)
    reportfile.write("run start time = " + time_string + "\n")
    frame_counter = D[first_sample:last_sample,1]
    mce_counter = D[first_sample:last_sample, 2]
    unix_time = (D[first_sample:last_sample, 3] - D[first_sample,3]) # not sure why the stragne scaling
    
    tmp= np.polyfit(T, unix_time,  1)
    smurf_to_unix_slope = tmp[0]
    print("slope = " + str(smurf_to_unix_slope))
    reportfile.write("unix time / smurf time slope = " + str(smurf_to_unix_slope)  + "\n")
    
    tmp = np.polyfit(T, frame_counter, 1)
    smurf_frame_rate = tmp[0]
    print("smurf frame rate = " + str(smurf_frame_rate ))
    reportfile.write("SMURF frame rate = " + str(smurf_frame_rate)  + "\n")
    
    output_rate = 1.0/np.mean(np.diff(T))
    print("average output rate = " + str(output_rate))
    reportfile.write("output rate = " + str(output_rate)  + "\n")
  
    tmp = np.polyfit(T, mce_counter, 1)
    mce_rate = tmp[0] 
    
    print("syncbox rate = " + str(mce_rate))
    reportfile.write("syncbox rate = " + str(mce_rate) + "\n")
    

    deltamce = np.diff(mce_counter)
    deltamce = np.append(deltamce, deltamce[-1])
    deltaframe = np.diff(frame_counter)
    deltaframe = np.append(deltaframe, deltaframe[-1])
    
    if R['plot_raw']:
        plt.figure(50)
        for n in range(0, num_pairs):  # loop to find elnods
            plt.plot(D[:,2*n+timing_channels] - np.mean(D[:,2*n+timing_channels]))
            plt.plot(D[:,2*n+timing_channels+1] - np.mean(D[:,2*n+timing_channels+1]))
        plt.title('raw data overlay')
        plt.show(block=False)
        
        
    Td = np.diff(T)  # difference milliseconds between outputs 
    Td = np.append(Td, Td[-1])
    fs = output_rate # average sample frequency

    wmax = int(np.floor(np.log((samples / welch_blocks)) / np.log(2)))
    welch_segments = np.power(2,wmax)    
    print("welch samples = ", welch_segments)
    f, Tx = signal.welch(Td, fs, nperseg = welch_segments)
    Tw = np.sqrt(Tx)
    print("sample frquency = ", fs)
      
    Rx = np.sign(1-np.mod(np.linspace(1, samples, samples)-.1 ,welch_segments))  # break up into timign steps for processing
    idx = np.where(Rx > 0)[0]
    if R['skip_ends']:
        Rx[idx[0]] = -1
        Rx[idx[-1]] = -1
    if R['show_timing_plots']:
        plt.figure(1000) # difference plot
        plt.plot(T, Td, '.')
        plt.xlabel('time seconds')
        plt.ylabel('delta time from Smurf timing')
        plt.title('time between outptut frames measured by SMRUF timing' + title)
        plt.show(block=False)
    
        plt.figure(1001)
        plt.loglog(f, Tw, 'k-')
        plt.title('Power spectrum of variation in delay between frames' + title)   
        plt.show(block=False)

        plt.figure(1002)
        plt.plot(T, deltaframe, 'k.')
        plt.title('Variaion of frames between outputs' + title)
        plt.show(block=False)   
    
        plt.figure(1003)
        plt.plot(T, deltamce, 'k.')
        plt.title('Variation in MCE frames between outputs' + title)
        plt.show(block=False)
 
        plt.figure(1004)
        plt.title('test plot')
        plt.plot(D[:, timing_channels])

        plt.figure(1005)
        plt.title('Rx diagnostic')
        plt.plot(Rx)

        
    if num_pairs > max_pairs:
        num_pairs = max_pairs
    
    er = np.zeros(num_pairs)    
    for p in range(0, num_pairs):
        A[:,p] = D[first_sample: last_sample, 2*p+timing_channels]
        B[:,p] = D[first_sample: last_sample, 2*p+timing_channels+1]
        Ma[p] = np.mean(A[:,p])
        Mb[p] = np.mean(B[:,p])
        Ax = A[:,p] - Ma[p]
        Bx = B[:,p] - Mb[p]
        Tx = D[first_sample:last_sample, 0] # timing
        er[p] = 1 # if no elnod
        
        
        f, Aw = signal.welch(Ax, fs, nperseg = welch_segments)
        f, Bw = signal.welch(Bx, fs, nperseg = welch_segments)
        
        (smurf_chanA, bandA, bchanA) = Cmap.get_channel(2*p)
        (smurf_chanB, bandB, bchanB) = Cmap.get_channel(2*p+1)
        gcpA = Cmap.gcp_index(2*p)
        gcpB = Cmap.gcp_index(2*p+1)
        chstrA = "A=" + str(smurf_chanA) + " b=" + str(bandA) + " ch=" + str(bchanA) 
        chstrB = "B=" + str(smurf_chanB) + " b=" + str(bandB) + " ch=" + str(bchanB)
        if R['use_gcp_names']:
            chstrA = chstrA + " gcp=" + str(gcpA)
            chstrB = chstrB + " gcp=" + str(gcpB)    
        title2 = R['run_name'] +' ' + R['pair_name'][p]+ '\n' + " pair = " + str(p) + " " 
        if enable_time_plot:


          
            fignum = p
            plt.figure(fignum)
            plt.clf()
            plt.plot(Tx, Ax, 'r-')
            plt.plot(Tx, Bx, 'b-')
            ec = 1
            if R['use_elnod']:
                ec= 0
                ne = len(R['elnod_starts'])
                for rs in R['elnod_starts']:
                    sp = rs + R['elnod_width']
                    #Ae = Ax[rs:sp]
                    #Be = Bx[rs:sp]
                    Ae = D[rs:sp,2*p+timing_channels]
                    Be =  D[rs:sp,2*p+timing_channels+1]
                    ec=  ec + elnodcalib(Ae, Be, R['elnod_lowpass']) / ne #normalized
            Cx = Ax - Bx * ec # elnod corrected
            if R['use_elnod']:
                plt.plot(Tx, Cx, 'g-')

            b,a = signal.butter(4, 1.0/R['lowpass_pts']) # butterworth filter generation
            Af = signal.filtfilt(b, a, Ax)
            Bf = signal.filtfilt(b, a, Bx)
            Cf = signal.filtfilt(b, a, Cx)
            if R['use_elnod']:
                plt.plot(Tx, Cf, 'c-')
            
            plt.grid()
            plt.xlabel(R['x_time_name'])
            plt.legend([chstrA, chstrB, "elnod diff", "diff filtered" ],  loc='lower left')
            plt.title(title2)
            if show_plots:
                plt.show(block=False)
            plt.savefig(R['plot_directory'] + run_string + R['plot_suffix'] +'_time_data_' + str(fignum), bbox_inches='tight')
        
        if R['show_welch_plots']:
            plt.figure(p+100)  # new figure
            plt.loglog(f, np.sqrt(Aw), 'r-')
            plt.loglog(f, np.sqrt(Bw), 'b-')
            plt.title(title2)
            plt.xlim(R['plot_fmin'], R['plot_fmax'])
            plt.xlabel(R['x_frequency_name'])
            plt.ylabel('raw counts/sqrt(Hz)')
            if R['show_plots']:
                plt.show(block=False)
            plt.savefig(R['plot_directory'] + run_string + R['plot_suffix']+ '_frequency_data_' + str(fignum), bbox_inches='tight')
        print(" p = " + str(p))
        # Next section plots correlation plots
        (freqs, apsx, pfit) = process_data(T, Ax, Bx, Rx, R)  # process data
        plot_data(freqs, apsx, pfit, p+2000, title2, run_string, R, R['aw_per_count'], "aW/sqrt(hz)", chstrA, chstrB)
        plot_data(freqs, apsx, pfit, p+3000, title2, run_string, R, R['pa_per_count'], "pA/sqrt(hz)", chstrA, chstrB)


        """
process-data takes an input data set and returns the correlated and uncorrelated data fourier transform
freqs is frequency values
aspx are the 4 arrasy,  A,  B,  A after subtracting fit to b,  B after subtracting fit to a
pfit is a commplex number of the scale of A / B per frequency. 
IMPORTANT!!!! Scaling applied to get corect single channle noise - but changes gain for physics signals (by sqrt(2))
IMPORTANT!!!!  See line with label "SCALING" 
DO NOT IGNORE ABOVE - THIS IS IMPORTANT!!!!

"""


def elnodcalib(A, B, lowpass): # do elnod calibration
    b,a = signal.butter(4, 1.0/lowpass) # butterworth filter generation
    Af = signal.filtfilt(b, a, A)
    Bf = signal.filtfilt(b, a, B)
    print('doing elnod')
    
    try:
        z = np.polyfit(Bf, Af, 1) # poly fit
    except:
        z = [1,1]
    return(z[0])
    
    
def process_data(tm, ax, bx, rx, R): # sx is the scale, defaults to 1
    ne = len(R['elnod_starts'])
    ec = 0
    if R['use_elnod']:
        ec= 0
        for rs in R['elnod_starts']:
            sp = rs + R['elnod_width']
            Ae = ax[rs:sp]
            Be = bx[rs:sp]
            ec=  ec + elnodcalib(Ae, Be, R['elnod_lowpass']) / ne #normalized

    cx = ax - bx * ec # elnod corrected
    pfit_num = R['poly_order']; # 3rd order polynomial fit for data
    idx = np.where(rx > 0)[0] # find markers - ugly index mess - i hate python
    nrec = idx.size # number of records
    record_pts = idx[1] - idx[0]  # for now assume equal length records (probaly OK, but ugly)
    h = np.hamming(record_pts)
    #h = signal.blackman(record_pts)
    hm = h / np.sqrt(np.mean(h*h))  # RMS normalized hamming window
    tm = tm[idx[0]:idx[1]]  # time 
    tstep = tm[1]-tm[0]  # time step
    freqs = np.fft.fftfreq(record_pts, tstep)  # list fo frequencies
    fpts = int(record_pts / 2)  # postitive frequecy points
    afft = np.zeros((record_pts, nrec,5), dtype=np.complex_) # holds fft data for 2 channels and 2 subracted channels
    pfit = np.zeros((record_pts),dtype=np.complex_)  #fit for slope
    for n in range(0,3):  # process both channels. (ugly!!)
        a = np.zeros((record_pts, nrec))     # create array for record
        af = np.zeros((record_pts, nrec))    # create array poly fit subtracted
        afh = np.zeros((record_pts, nrec))   # create array for hamming window multilied
        if n == 0:              #THIS IS UGLY UGLY, too lazy to do correctly. 
            x = ax
        elif n ==1:
            x = bx
        else:
            x = cx
        for j in range(0,nrec-1):               # loop over records for spectrum
            a[:,j] = x[idx[j]:idx[j+1]]   
            P = np.polyfit(tm, a[:,j],pfit_num) # polyfit to remove trends
            r = np.polyval(P,tm)
            af[:,j] = a[:,j]-r
            afh[:,j] = af[:,j]* hm              #multiply by hamming windows
            # mystery 2.5 factor that I dno't understand. very bad
            afft[:,j,n] = np.fft.fft(afh[:,j]) * np.sqrt(tstep) / np.sqrt(fpts) # do fft (note complex out)
    for r in range(0, fpts):  # fit by frequency band
        at = afft[r,:,0]
        bt = afft[r,:,1]
        a_s = np.sum(np.power(np.abs(at),2))  #power of A channel in this band
        b_s = np.sum(np.power(np.abs(bt),2))
        try:
            P = np.polyfit(bt, at, 2)  #poly fit a chnnel to b channel in preperation for subtraction
            pfit[r] = P[1]
            cta = np.polyval(P, bt) # Prediction of A based on B
            c_sa = np.sum(np.power(np.abs(cta),2)) # Power of fit data in this band
            va = np.sqrt(1+c_sa / (a_s +c_sa)) # SCALING of spectrum to correct based on how much correlated data is included
            dt = at-cta # subtracted fit
            afft[r,:,3] = dt / va # copy into difference fft
                                                # Same as above for B channel
            P = np.polyfit(at, bt, 2)
            ctb = np.polyval(P, at) 
            c_sb = np.sum(np.power(np.abs(ctb),2)) # rms sums
            vb = np.sqrt(1+c_sb / (b_s +c_sb))
            dt = bt-ctb # subtracted fit
            afft[r,:,4] = dt / vb # copy into difference fft
        except:  # if on array is all zeros or nan or somethikg
            afft[r,:,3] = 0 * afft[r,:,0]
            afft[r,:,3] = 0 * afft[r,:,0]
    apsx = np.zeros((record_pts, 5)) # holds output data
    for n in range(0,5):
        for j in range(0, nrec-1):
            apsx[:,n] = apsx[:,n]  + np.power(np.abs(afft[:,j,n]),2) / (nrec+1)
    return(freqs[0:fpts], apsx[0:fpts,:], pfit[0:fpts])


def plot_data(freqs, apsx, pfit, fignum,title,  fname_prefix, R, scl, vtxt, chstrA, chstrB):
     plt.figure(fignum)
     plt.clf()

     plt.loglog(freqs, scl * np.sqrt(apsx[:,0]), 'r-')  # channel a
     plt.loglog(freqs, scl * np.sqrt(apsx[:,1]), 'b-')  # chanmnel b
              
     if R['show_decorrelated']:
         plt.loglog(freqs, scl * np.sqrt(apsx[:,3]), 'm-')  # a after removing contribution from b
     if R['use_elnod']:
         plt.loglog(freqs, (1.0 /np.sqrt(2)) * scl * np.sqrt(apsx[:,2]), 'g-')   # simple subtraction
     leg0 = [chstrA, chstrB]
     if R['show_decorrelated']:
         leg0.append("Chan A uncorrelated")
     if R['use_elnod']:
         leg0.append("Elnod corrected")
    
     plt.legend(leg0)
     plt.xlabel("Frequency Hz")
     plt.ylabel(vtxt)
    
     plt.xlim(R['plot_fmin'], R['plot_fmax'])
     plt.grid( which = 'major')
     plt.grid( which = 'minor')
     plt.title(title)    
     if R['show_plots']:
         plt.show(block=False)
     plt.savefig(R['plot_directory'] + R['run_name'] + R['plot_suffix']+ '_'+'corr_data_' + str(fignum), bbox_inches='tight') # disable for now

