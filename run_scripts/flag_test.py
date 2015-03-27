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
Z_ant = sqrt(R_ant**2+X_ant**2)
Z_ant_sm = itp.UnivariateSpline(F_ant,Z_ant,s=0.01)
Effic = ef.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic,s=0.01)
Eff_50 = ef.effic(50.*ones(len(F_amp)),zeros(len(F_amp)),F_amp,R_amp,X_amp,F_amp)
Eff_50_sm = itp.UnivariateSpline(F_amp,Eff_50,s=0.01)
Eff_100 = ef.effic(100.*ones(len(F_amp)),zeros(len(F_amp)),F_amp,R_amp,X_amp,F_amp)
Eff_100_sm = itp.UnivariateSpline(F_amp,Eff_100,s=0.01)

total_mask = 0.
mid_mask = 0.
spike_m = 0.
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
    short_data = loadtxt(Kdir+'June_'+date_ind+'_fit_short.txt')
    load_data = loadtxt(Kdir+'June_'+date_ind+'_fit_50ohm.txt')
    term_data = loadtxt(Kdir+'June_'+date_ind+'_fit_100ohm.txt')
    noise_data = loadtxt(Kdir+'June_'+date_ind+'_fit_noise.txt')

    short_full = loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')
    load_full = loadtxt(Kdir+'June_'+date_ind+'_avg_50ohm.txt')
    term_full = loadtxt(Kdir+'June_'+date_ind+'_avg_100ohm.txt')
    noise_full = loadtxt(Kdir+'June_'+date_ind+'_avg_noise.txt')


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
                int_mask = mask
                int_mask = ff.flagging(sub_data,freq,2.,freqscale)
                mid_mask = mid_mask+sum(int_mask)
                total_sum = total_sum +len(int_mask)

#adding in freq spike flagging
                mask = ff.spike_flag(sub_data,int_mask,freq,15.)
                spike_m = spike_m+sum(mask)

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
print 'Percent of Data Flagged with Spike Flagger:',100.*spike_m/total_sum

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


#Adding in the short flagging:

for i in range(0,len(sortdata)):
    new_mask = ff.cal_flag(short_full,short_data,sortmask[i],freq,1.e-10)
    sortmask[i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data adding short masking',percent_masked


#Adding in the threshold flagging:
for i in range(0,len(sortdata)):
    new_mask = ff.threshold_flag(sortdata[i],sortmask[i],freq,5.)
    sortmask[i] = new_mask 
 
percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0])) 
print 'Percentage of Masked Data adding threshold masking',percent_masked 


fmin = where(freq<=50.)[0][-1]
fmax = where(freq<=100.)[0][-1]

percent_masked_trunc = 100.*sum(sortmask[:,fmin:fmax])/(len(sortmask)*(fmax-fmin))
print 'Percentage of Masked Data for ',freq[fmin],' to ',freq[fmax],' MHz: ',percent_masked_trunc


#Calculating time mean
mean_sd, mean_sm = cf.time_mean(sortdata,sortmask)

Kt_mask = ff.spike_flag(mean_sd,mean_sm,freq,10.)
print 'Percent of mean data masked',100.*sum(Kt_mask)/len(Kt_mask)

#Applying calibration data
Kt = 300./(load_data-short_data)
mean_Kt = Kt*(mean_sd-short_data)*Eff_sm(freq) 

KA = (load_data-term_data)/(Eff_50_sm(freq)-4*Eff_100_sm(freq))
Kt50 = 300./(load_data-short_data-KA*Eff_50_sm(freq))
Kt100 = 300./(term_data-short_data-KA*4*Eff_100_sm(freq))
PZ = KA*(Z_ant_sm(freq)**2)*Eff_sm(freq)/50**2
mean_Kt50 = Kt50*(mean_sd - short_data-PZ)*Eff_sm(freq)

pylab.scatter(freq, mean_Kt,c='b',edgecolor='b',s=3)
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt,mask=Kt_mask)),'c--')
pylab.scatter(freq, mean_Kt50, c= 'r',edgecolor='r',s=3)
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt50,mask=Kt_mask)),'m--') 
pylab
pylab.xlim(50,120)
pylab.xlabel('Frequency (MHz)')
pylab.ylim(0,6e3)
pylab.ylabel('Temperature (Kelvin)')
pylab.grid()
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_comp',dpi=300)
pylab.clf()


for i in range(0,len(sortdata)):
    for j in range(0,len(sortdata[0])):
        if sortmask[i,j] == 1.0:
            sortdata[i,j] = 0.0
        else:
            sortdata[i,j] = Kt50[j]*(sortdata[i,j] - short_data[j]- PZ[j])*Eff_sm(freq[j])

pylab.imshow(sortdata[:,fmin:fmax],vmin=0,vmax=6e3,aspect=50./(sorttime[-1]-sorttime[0]),extent=(50,100,sorttime[-1],sorttime[0]))
pylab.title('Masked Data')
pylab.colorbar()
pylab.title('Measured Sky Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.savefig(outdir+'June_'+date_ind+'_update_jncal_masked_waterfall',dpi=300)
pylab.clf()
