"""
Module to remove foreground signal using polynomial fitting
for a single day of data. 
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
import numpy.polynomial.polynomial as poly
import file_funcs as ff
import eff_funcs as ef
import cal_funcs as cf

#Main directories for the input and output
indir = '/lustre/tcv/freq_rebinned_data/'
outdir = '/lustre/tcv/mean_cal_data/'
Kdir = '/lustre/tcv/calibration_data/'
directories = os.listdir(indir)

#Setting a single day for parallel computation
date_ind = sys.argv[1]
minf = 60.
maxf = 90.

#Additional files needed
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

# Efficiency calculation
R_ant,X_ant,F_ant = ef.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = ef.imped_skrf(amp_s_file,0.0)
Effic = ef.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic,0.01)

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
 
#Prepping calibration data for use
    short_data = loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')
    load_data = loadtxt(Kdir+'June_'+date_ind+'_avg_50ohm.txt')
    K_gsm = loadtxt(Kdir+'June_'+date_ind+'_K_gsm.txt')
    K_dgsm = loadtxt(Kdir+'June_'+date_ind+'_K_dgsm.txt')

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

#Calculating time mean
mean_sd, mean_sm = cf.time_mean(sortdata,sortmask)

#Additional Frequency Masking for Mean Data
percent_masked = 100.*sum(mean_sm)/(len(mean_sm))
print 'Percent of mean data masked', percent_masked

spike_flag = ff.spike_flag(mean_sd,mean_sm,freq,5.)
for i in range(0,len(mean_sd)):
    if mean_sm[i]==1.0:
        print 'Mean Mask Frequency at:',freq[i]
    if spike_flag[i]==1.0:
        mean_sm = 1.0
        print 'Spike Frequency at:', freq[i]    

#Applying calibration data
Kt = 300./(load_data-short_data)

mean_Kt = Kt*mean_sd
mean_Kdgsm = K_dgsm*mean_sd
mean_Kgsm = K_gsm*mean_sd

#Foreground fitting - only going to do n=2 for now.
Kt_fit, Kt_params = cf.poly_fore(mean_Kt,mean_sm,freq,minf,maxf,2)
Kdgsm_fit,Kdgsm_params = cf.poly_fore(mean_Kdgsm,mean_sm,freq,minf,maxf,2)
Kgsm_fit,Kgsm_params = cf.poly_fore(mean_Kgsm,mean_sm,freq,minf,maxf,2)

#Foreground parameter comparison:
print 'Foreground parameters for Kt are:',Kt_params
print 'Foreground parameters for Kgsm are:',Kgsm_params
print 'Foreground parameters for Kdgsm are:',Kdgsm_params

#Outputting fit and mean data
savetxt(outdir+'June_'+date_ind+'_mean_mask.txt',mean_sm,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_Kt_mean.txt',mean_Kt,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_Kt_fit2.txt',Kt_fit,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_Kgsm_mean.txt',mean_Kgsm,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_Kgsm_fit2.txt',Kgsm_fit,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_Kdgsm_mean.txt',mean_Kdgsm,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_Kdgsm_fit2.txt',Kdgsm_fit,delimiter=' ')

#Comparison plot
pylab.scatter(freq,mean_Kt,color='b',edgecolor='b',label='Kt mean')
pylab.plot(freq,Kt_fit,c='b',label='Kt fit')
pylab.scatter(freq,mean_Kgsm,color='g',edgecolor='g',label='Kgsm mean')
pylab.plot(freq,Kgsm_fit,c='g',label='Kgsm fit')
pylab.scatter(freq,mean_Kdgsm,color='r',edgecolor='r',label='Kdgsm mean')
pylab.plot(freq,Kdgsm_fit,c='r',label='Kdgsm fit')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.legend()
pylab.grid()
pylab.title('Polynomial Fit Comparisons')
pylab.savefig(outdir+'June_'+date_ind+'_mean_fits',dpi=300)
pylab.clf()
