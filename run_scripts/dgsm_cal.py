"""
Module to calculate the K_dgsm for a single day of data.
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
import ephem as eph
import matplotlib.pyplot as plt
import file_funcs as ff
import eff_funcs as ef
import cal_funcs as cf

#Main directories for the cal and input data
caldir = '/lustre/tcv/truncated_data/'
indir = '/lustre/tcv/freq_rebinned_data/'
Kdir = '/lustre/tcv/calibration_data/'
directories = os.listdir(maindir)

#Setting a single day for parallel computation
date_ind = sys.argv[1]

#Additional files needed
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_update.dat')


#Setting rebinning scales
timescale = 32
freqscale = 32

# Efficiency calculation
R_ant,X_ant,F_ant = ef.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = ef.imped_skrf(amp_s_file,0.0)
Effic = ef.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic,0.01)

#GSM data rearranged to array with first index = time, second index = freq
gsm_freq = arange(60,90,1)
gsm_time = arange(0,24,24./289)
gsm_data = zeros((len(gsm_time),len(gsm_freq)))
for t in range(0,len(gsm_time)):
    for f in range(0,len(gsm_freq)):
        gsm_data[t,f] = gsm_raw_data[f*289+t,2]

# Make data array
processed_data = []
processed_mask = []
processed_time = []
processed_mtime = []

if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being analyzed is:',direct
    directory = maindir+direct+'/'
    dirlist = os.listdir(directory)
 
#Prepping short data for use
    short_data = loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')

#Iterate for each file in the directory
   for fname in dirlist:
        if len(fname.split('_'))>=3:
            filename = directory+fname
#load data file
            time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40.+width/2,130.,width)
            if len(freq)>len(sub_data):
                freq = freq[0:-2]
            if fname.split('_')[-1]=='mask.dat':
                if fname.split('_')[-2]=='antenna':
                    processed_mask.append(sub_data)
                    processed_mtime.append(time)
            elif fname.split('_')[-1]=='antenna.dat':
                new_data = (sub_data/Eff_sm(freq)-short_data)
                processed_data.append(new_data)
                processed_time.append(time)

#Setting up sidereal time arrays
sid_time = []
sid_mtime = []

idate = '2013/6/1'
for i in range(0,len(processed_time)):
    single_sid = ff.sidereal(processed_time[i],idate)
    sid_time.append(single_sid)
for i in range(0,len(processed_mtime)):
    single_sid = ff.sidereal(processed_time[i],idate)
    sid_mtime.append(single_sid)

#Sort data by LST
processed_data = array(processed_data)
processed_mask = array(processed_mask)

sortind = argsort(sid_time)
sorttime = zeros(len(processed_data))
sortdata = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(sid_time)):
    sorttime[i] = sid_time[sortind[i]]
    sortdata[i] = processed_data[sortind[i]]
sortindm = argsort(sid_mtime)
sortmask = zeros((len(processed_data),len(processed_data[0])))
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

#Time rebinning to match gsm binning (5 minutes)
stack_data, stack_mask = cf.match_binning(gsm_time,freq,sorttime,sortdata,sortmask)

#Only grabbing the gsm times that we have measured data for:
lim_stack,lim_mask,lim_gsm,lim_time = cf.lim_bin(freq,stack_data,stack_mask,gsm_freq,gsm_data,gsm_time)

#Calculating time means of both GSM and Data
mean_sd, mean_sm = cf.time_mean(lim_stack,lim_mask)
gsm_mask = zeros(len(lim_gsm))
mean_gsm, mean_gmask = cf.time_mean(lim_gsm,gsm_mask)

#Generating mean subtracted arrays
lim_ms_stack= zeros((len(lim_stack),len(stack_data[0])))
for i in range(0,len(lim_stack)):
    lim_ms_stack[i] = lim_stack[i]-mean_sd
lim_ms_gsm = zeros((len(lim_stack),len(freq)))
for i in range(0,len(lim_stack)):
    lim_ms_gsm[i] = lim_gsm[i]-mean_gsm

lim_ms_stack = array(lim_ms_stack)
lim_ms_gsm = array(lim_ms_gsm)

#Calculating Gain (K_dgsm)
K_dgsm = []
K_dgsm0 = [1e10]
for f in range(0,len(freq)):
    if mean_mask[f]==1.0:
        K_dgsm.append(K_dgsm0)
        print 'Fitting failed at frequency', freq[f]
    else:
        Kf = cf.gain_calc(lim_ms_stack[:,f],lim_mask[:,f],lim_ms_gsm[:,f],K_dgsm0)
        K_dgsm0 = Kf
        K_dgsm.append(Kf)
    
K_dgsm = array(K_dgsm)
K_dgsm = abs(K_dgsm)
f70 = where(freq<=70.)[0][-1]
print 'K_dgsm at 70 MHz is:',K_dgsm[f70]

savetxt(Kdir+'June_'+date_ind+'_K_dgsm.txt',K_dgsm,delimiter=' ')

#Plot results
pylab.plot(freq,mean_gsm,label='mean GSM Spectrum')
pylab.plot(freq,mean_sd*K_dgsm,label='mean Calibrated Spectrum')
pylab.grid()
pylab.xlim(40,130)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Time Averaged Signals')
pylab.legend()
pylab.savefig(Kdir+'June_'+date_ind+'_K_dgsm_mean_spectrum',dpi=300)
pylab.clf()

pylab.plot(lim_time,lim_ms_gsm[:,f70],label='GSM time series')
pylab.plot(lim_time,lim_ms_stack[;,f70]*K_dgsm[f70],label='Calibrated time series')
pylab.xlim(0,24)
pylab.xlabel('Local Sidereal Time (Hours)')
pylab.ylabel('Temperature (Kelvin)')
pylab.grid()
pylab.title('Time Dependence at 70 MHz')
pylab.legend()
pylab.savefig(Kdir+'June_'+date_ind+'_K_dgsm_time_series',dpi=300)
pylab.clf()
