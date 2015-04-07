"""
Module to compile data for plotting.
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
import os
import skrf as rf
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('../../hibiscus'))
import file_funcs as ff

#Main directories for the input and output
#indir = '/lustre/tcv/truncated_data/'
#indir = '../../Reverb_chamber_data/'
#indir = '../../Karoo_data/truncated/'
indir = '../../Karoo_data_day1/rebinned/'
#outdir = '/lustre/tcv/freq_rebinned_data/'
outdir = '../../Karoo_data_day1/rebinned/'
#Kdir = '/lustre/tcv/rfi_check_data/'
#directories = os.listdir(indir)

#Setting a single day for parallel computation
#date_ind = sys.argv[1]

#Setting Rebinning scales
timescale = 32
freqscale = 32

processed_data = []
processed_mask = []
processed_time = []
processed_mtime = []


new_directory = outdir
#if freqscale==32:
#repping calibration data for use
#    short_data = loadtxt(Kdir+'June_'+date_ind+'_fit_short.txt')
#    short_full = loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')

#    dirlist = os.listdir(directory)
dirlist = os.listdir(indir)
    #print dirlist
#Iterate for each file in the directory
for fname in dirlist:
    if len(fname.split('-'))>=5:
        if fname.split('-')[0]=='2015':
#            if fname.split('-')[1]=='.dat':
#            if fname.split('-')[-1]!='cal.dat':
#                filename = directory+fname
            filename = indir + fname
#            if os.path.isfile(indir+'rebinned/'+fname)==False:
#                print fname
#load data file
            time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40,130.0,width)
            if len(freq)>368:
                print 'problem', fname
            elif len(freq)<367:
                print 'problem',fname
            
            if fname.split('_')[-1]=='mask.dat':
                if fname.split('_')[-2]=='antenna':
                    processed_mask.append(sub_data)
                    processed_mtime.append(time)
            elif fname.split('_')[-1]=='antenna.dat':
                processed_data.append(sub_data)
                processed_time.append(time)

print shape(processed_data)
print shape(processed_mask)
processed_time = processed_time - 9.*ones(len(processed_time))
processed_mtime = processed_mtime - 9.*ones(len(processed_time))

sid_time = []
sid_mtime = []

idate = '2015/4/1'
longitude = '21.4109'
latitude = '-30.7216'
for i in range(0,len(processed_time)):
    single_sid = ff.sidereal(processed_time[i],idate,longitude,latitude)
    sid_time.append(single_sid)
for i in range(0,len(processed_mtime)):
    single_sid = ff.sidereal(processed_mtime[i],idate,longitude,latitude)
    sid_mtime.append(single_sid)

processed_data = array(processed_data)
processed_mask = array(processed_mask)

sortind = argsort(processed_time)
sorttime = zeros(len(processed_data))
sortdata = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(sid_time)):
    sorttime[i] = processed_time[sortind[i]]
    sortdata[i] = processed_data[sortind[i]]
sortindm = argsort(processed_mtime)
sortmask = zeros((len(processed_mask),len(processed_mask[0])))
for i in range(0,len(sid_mtime)):
    sortmask[i] = processed_mask[sortindm[i]]

#Frequency Masking Amount Calculation
percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency Masking',percent_masked

#Adding in Time Masking
for i in range(0,len(freq)):
    new_mask = ff.timeflag(sortdata[:,i],sortmask[:,i],sorttime,3.,timescale)
    sortmask[:,i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking',percent_masked

