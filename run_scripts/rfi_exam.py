"""
Module to rebin a single day of data
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
sys.path.append(os.path.abspath('/home/tcv/hibiscus'))
import file_funcs as ff
import eff_funcs as ef
import cal_funcs as cf

#Main directories for the input and output
indir = '/lustre/tcv/truncated_data/'
outdir = '/lustre/tcv/rfi_check_data/'
Kdir = outdir
directories = os.listdir(indir)

#Setting a single day for parallel computation
date_ind = sys.argv[1]

#Setting Rebinning scales
timescale = 32
freqscale = 32

ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'
# Efficiency calculation
R_ant,X_ant,F_ant = ef.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = ef.imped_skrf(amp_s_file,0.0)
Effic = ef.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic,s=0.01)


total_mask = 0.
mid_mask = 0.
total_sum = 0.

processed_data = []
processed_mask = []
processed_time = []
processed_mtime = []

if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being rebinned is:',direct
    directory = indir+direct+'/'
    new_directory = outdir+direct+'/'
    dirlist = os.listdir(directory)

#repping calibration data for use
    short_data = loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')
    load_data = loadtxt(Kdir+'June_'+date_ind+'_avg_50ohm.txt')

#Iterate for each file in the directory
    for fname in dirlist:
        if len(fname.split('-'))>=3:
            if fname.split('_')[-1]=='antenna.dat':
                filename = directory+fname
#load data file
                time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 90.0/len(sub_data)
                freq = arange(40,130.0,width)
#basic freq flagging
                mask = ff.flagging(sub_data,freq,3.,freqscale)
                mid_mask = mid_mask+sum(mask)
                total_sum = total_sum +len(mask)

#Check for nan/inf (should be nulls)
                nandata = where(isnan(sub_data))[0]
                for i in range(0,len(nandata)):
                    sub_data[nandata[i]] = 0.0
                infdata = where(isinf(sub_data))[0]
                for i in range(0,len(infdata)):
                    sub_data[infdata[i]] = 0.0

                processed_data.append(sub_data)
                processed_mask.append(mask)
                processed_time.append(time)
                processed_mtime.append(time)

print 'Percent of Data Flagged without Spike Flagger:',100.*mid_mask/total_sum

#Setting up sidereal time arrays
sid_time = []
sid_mtime = []

idate = '2013/6/1'
for i in range(0,len(processed_time)):
    single_sid = ff.sidereal(processed_time[i],idate)
    sid_time.append(single_sid)
for i in range(0,len(processed_mtime)):
    single_sid = ff.sidereal(processed_mtime[i],idate)
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

#Calculating time mean
mean_sd, mean_sm = cf.time_mean(sortdata,sortmask)

#Additional Frequency Masking for Mean Data
percent_masked = 100.*sum(mean_sm)/(len(mean_sm))
print 'Percent of mean data masked', percent_masked

#Applying calibration data
Kt = 300./(load_data-short_data)
mean_Kt = Kt*mean_sd

Kt_mask = ff.spike_flag(mean_Kt,mean_sm,freq,10.)
print 'Percent of Kt mean data masked',100.*sum(Kt_mask)/len(Kt_mask)

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_sd,mask=Kt_mask)))
pylab.xlim(60,120)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Data before calibration')
pylab.ylim(0,6e-8)
pylab.grid()
pylab.title('Time Mean without calibration')
#pylab.savefig(outdir+'June_'+date_ind+'_mean_uncal_spectrum',dpi=300)
pylab.clf()

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt,mask=Kt_mask)))
pylab.xlim(60,120)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(0,6e3)
pylab.grid()
pylab.title('Time Mean with JNC calibration')
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_spectrum',dpi=300)
pylab.clf()

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt,mask=Kt_mask)))
pylab.xlim(80,100)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(1e3,3e3)
pylab.grid()
pylab.title('Time Mean with JNC calibration')
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_spectrum_zoom',dpi=300)
pylab.clf()


fmin = where(freq<=60.)[0][-1]
fmax = where(freq<=120.)[0][-1]
for i in range(0,len(sortind)):
    sortdata[i] = sortdata[i]*Kt

pylab.imshow(sortdata[:,fmin:fmax],vmin=0,vmax=6e3,aspect=60./len(sortind),extent=(60,120,len(sortind),0.0))
pylab.colorbar()
pylab.title('Variation of data over the day (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (sample index)')
pylab.savefig(outdir+'June_'+date_ind+'_unmasked_cal_waterfall',dpi=300)
pylab.clf()

for i in range(0,len(sortind)):
    for j in range(0,len(freq)):
        if sortmask[i,j]==1:
            sortdata[i,j] = 0

pylab.imshow(sortdata[:,fmin:fmax],vmin=0,vmax=6e3,aspect=60./len(sortind),extent=(60,120,len(sortind),0.0))
pylab.colorbar() 
pylab.title('Variation of data over the day (Kelvin), masked data == zeros')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (sample index)')
pylab.savefig(outdir+'June_'+date_ind+'_masked_cal_waterfall',dpi=300)
pylab.clf()

