# -*- coding: utf-8 -*-
"""
procdata
Script to call analysis routines 

"""
import subprocess
import analyzedata
import file_reader

#defaultdatadir = '/home/frisch/User4/Data/' # list here at top since it changes oftern
defaultdatadir = "/mnt/d/smurf_data/Pole_LMS_2/"


R0 = dict() # will hold base parameters to be updateed
R0['data_directory'] = defaultdatadir  #where the data is located
R0['temp_directory'] = '/tmp/'  # where the temporary text file is created
R0['plot_directory'] = './Plots/'  # where the plots and output file are saved
R0['run_name'] =  '' # prefix for all data
R0['description'] = "none" # short text header describing run
R0['list_file'] = 'pole_list.txt'   #default list of which channels in teh data file are read out
R0['mask_file'] = 'pole_mask.txt'  # default mask file
R0['report_file'] = 'report_file.txt'
R0['compile_command']= 'gcc analyzedata.cpp -o analyzedata'   
R0['extract_program'] = './analyzedata'  # the executible  compiled form the c code to extract text files
R0['num_averages'] = 0  #  this can average and downsample data
R0['chop_front'] = 0.0 # fraction of front of pulse to remove
R0['chop_back'] = 0.0 # fractino of back of pulse to remove
R0['welch_blocks'] = 20 # number of welch segements
R0['pa_per_count'] = 138.0  # just for now
R0['aw_per_count'] = 0.2*R0['pa_per_count']  # just for now.
R0['x_time_name'] = 'seconds '
R0['x_frequency_name'] = 'Hz '
R0['use_gcp_names'] = True
R0['txtfname'] = '/dev/null'  # just null for now, will be replaced 
R0['diag_timing'] = True   # use diagnostic timing outputs 
R0['cat_files'] = True   # concatinate all the standard smurf files into a single dat file.
R0['extract_txt'] = True # extracts to text file
R0['show_decorrelated'] = True
R0['use_elnod'] = False
R0['elnod_starts'] = []  # list of index for start of each elnoe.  (have to find this by hand)
R0['elnod_width'] = 0  # width of elnods - assumed the same for all
R0['elnod_lowpass'] = 100 # 100 points averaged for elnod - may help fit accuracy (or maybe not)
R0['lowpass_pts'] = 1000 # for low pass filter to display elnod corrected data
R0['plot_raw'] = False  # turn on to plot raw data - useful to find elnod points
R0['poly_order'] = 3  # for offset subtraction
R0['skip_ends'] = False # skip frist and last data segments
R0['show_timing_plots'] = False # dislay of various timing info, not usualy needed
R0['show_welch_plots'] = False # used as a diagnostics
R0['show_plots'] = True # actually display plots, also pauses between runs
R0['save_plots'] = True # save pngs of data
R0['plot_suffix'] = ''  # appended to end of plot names. 
R0['plot_fmin'] = 0.1 # minimum plot frequency Hz
R0['plot_fmax'] = 10 # maximum plot frequency Hz
R0['active'] = False  # set to true to actualy run this analysis
R0['fake_data'] = 0;  # fake data mode.    0 is don't fake data. 
R0['pair_name'] = []  # create an empty list
for n in range(0,2000):
    R0['pair_name'].append('')  # just append blank for now

Runlist = []  # creates a new list of jobs to run, only use more than one at a time if graphics is not displayed



""" THIS IS WHERE YOU START EDITING STUFF """
R = R0.copy()  # created an editable copy
R['pair_name'] = R0['pair_name'].copy()  # kludged deepcopy 
R['run_name']= '1569482492'  
R['chop_front'] = 0.05 # fraction of front of pulse to remove
R['chop_back'] = 0.05 # fractino of back of pulse to remove
R['welch_blocks'] = 100 # for long run
R['plot_fmin'] = .01 
R['description'] = 'Stanford '
R['pair_name'][0] = 'Transition'
R['pair_name'][1] = 'Superconducting'
R['list_file'] = 'stanfordtest_list.txt'
R['mask_file'] = 'stanfordtest_mask.txt'
R['num_averages'] = 0
R['use_gcp_names'] = True
R['cat_files'] = True  # concatinate all the standard smurf files into a single dat file.
R['extract_txt'] = True # extracts to text file 
R['active'] = False   # actually run this data set
Runlist.append(R)

R = R0.copy()  # created an editable copy
R['pair_name'] = R0['pair_name'].copy()  # kludged deepcopy 
R['run_name']= '1558717957'
R['description'] = 'Pole good running '
#R['chop_front'] = 0.05 # fraction of front of pulse to remove
#R['chop_back'] = 0.05 # fractino of back of pulse to remove
R['plot_fmin'] = 0.01 # minimum plot frequency
R['poly_order'] = 2  # for offset subtraction
R['welch_blocks'] = 5 # number of welch segements
R['use_elnod'] = True
R['elnod_starts'] = [3000, 489000]
R['elnod_width'] = 5500
R['skip_ends'] = True 
R['show_decorrelated'] = False
R['show_timing_plots'] = True
R['cat_files'] = False
R['active'] = False
Runlist.append(R)


R = R0.copy()  # created an editable copy
R['pair_name'] = R0['pair_name'].copy()  # kludged deepcopy 
R['run_name']= '1569280219'
R['description'] = 'Pole Data '
#R['chop_front'] = 0.05 # fraction of front of pulse to remove
#R['chop_back'] = 0.05 # fractino of back of pulse to remove
R['plot_fmax'] = 10 # maximum plot frequency
R['welch_blocks'] = 20 # number of welch segements
R['cat_files'] = False
R['active'] = False
Runlist.append(R)


R = R0.copy()  # created an editable copy
R['pair_name'] = R0['pair_name'].copy()  # kludged deepcopy 
#R['data_directory'] = "/mnt/d/smurf_data/Pole_LMS_2/"  # data directory on D
#R['data_directory'] = "/mnt/c/Data/"
R['run_name'] = '1571921948'
R['description'] = "Pole, LMS_gain = 2 "
R['plot_raw'] = True
R['elnod_starts'] = [1175000]
R['elnod_width'] = 175000
R['use_elnod'] = True
R['mask_file'] = R['data_directory'] + R['run_name'] + '_mask.txt'
R['plot_fmin'] = 0.1 # minimum plot frequency
R['plot_fmax'] = 100 # maximum plot frequency
R['plot_suffix'] = '_use_elnod_'
R['show_plots'] = True
R['active'] = True
Runlist.append(R)

R = R0.copy()  # created an editable copy
R['pair_name'] = R0['pair_name'].copy()  # kludged deepcopy 
#R['data_directory'] = "/mnt/d/smurf_data/Pole_LMS_2/"  # data directory on D
#R['data_directory'] = "/mnt/c/Data/"
R['run_name'] = '1571921948'
R['description'] = "Pole, LMS_gain = 2 full "
R['plot_raw'] = True
R['mask_file'] = R['data_directory'] + R['run_name'] + '_mask.txt'
R['chop_front'] = 0.0 # fraction of front of pulse to remove
R['chop_back'] = 0.0 # fractino of back of pulse to remove
R['plot_suffix'] = '_full_'
R['show_plots'] = False
R['active'] = False 
Runlist.append(R)

R = R0.copy()  # created an editable copy
R['pair_name'] = R0['pair_name'].copy()  # kludged deepcopy 
R['data_directory'] = "/mnt/d/smurf_data/pole_4kHz_unfiltered_vs_bias_20191009/"  # data directory on D
R['run_name'] = '1570622134'
R['description'] = "Pole, 4KHz "
R['plot_raw'] = True
R['mask_file'] = R['data_directory'] + R['run_name'] + '_mask.txt'
R['chop_front'] = 0.0 # fraction of front of pulse to remove
R['chop_back'] = 0.0 # fractino of back of pulse to remove
R['plot_fmin'] = 0.1 # minimum plot frequency
R['plot_fmax'] = 2000 # maximum plot frequency
R['active'] = False
Runlist.append(R)


""" THIS IS WHERE YOU STOP EDITING STUFF""" 



for R in Runlist:
    if not R['active']:
        continue
    print(R['run_name'])
    D = file_reader.file_reader(R)  # cats, and reads data file
    analyzedata.run(D, R)  # this does all the real work
    if R['show_plots']:
        input("Wait for plots to appear, then press enter to continue to the next run")
    

