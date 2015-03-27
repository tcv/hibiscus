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

print 'Percent of Data Flagged without Spike Flagger: ',100.*mid_mask/total_sum
print 'Percent of Data Flagged with Spike Flagger: ',100.*spike_m/total_sum

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
print 'Percentage of Masked Data from Frequency Masking: ',percent_masked

#Adding in Time Masking
for i in range(0,len(freq)):
    new_mask = ff.timeflag(sortdata[:,i],sortmask[:,i],sorttime,3.,timescale)
    sortmask[:,i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking: ',percent_masked

#Adding in the short flagging:

for i in range(0,len(sortdata)):
    new_mask = ff.cal_flag(short_full,short_data,sortmask[i],freq,1.e-10)
    sortmask[i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data adding short masking: ',percent_masked


#Adding in the threshold flagging:
for i in range(0,len(sortdata)):
    new_mask = ff.threshold_flag(sortdata[i],sortmask[i],freq,5.)
    sortmask[i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data adding threshold masking: ',percent_masked

fmin = where(freq<=50.)[0][-1]
fmax = where(freq<=100.)[0][-1]

percent_masked_trunc = 100.*sum(sortmask[:,fmin:fmax])/(len(sortmask)*(fmax-fmin))
print 'Percentage of Masked Data for ',freq[fmin],' to ',freq[fmax],' MHz: ',percent_masked_trunc

#Calculating time mean
mean_sd, mean_sm = cf.time_mean(sortdata,sortmask)

#Additional Frequency Masking for Mean Data
percent_masked = 100.*sum(mean_sm)/(len(mean_sm))
print 'Percent of mean data transfer masked: ', percent_masked


Kt_mask = ff.spike_flag(mean_sd,mean_sm,freq,10.)
Kt_mask = ff.cal_flag(short_full,short_data,Kt_mask,freq,1.e-10)
print 'Percent of mean data spike+cal masked: ',100.*sum(Kt_mask)/len(Kt_mask)



#Applying calibration data
Kt = 300./(load_data-short_data)
mean_Kt = Kt*(mean_sd-short_data)*Eff_sm(freq) 

KA = (load_data-term_data)/(Eff_50_sm(freq)-4*Eff_100_sm(freq))
Kt50 = 300./(load_data-short_data-KA*Eff_50_sm(freq))
Kt100 = 300./(term_data-short_data-KA*4*Eff_100_sm(freq))

#pylab.plot(freq,Eff_sm(freq))
#pylab.plot(freq,Eff_50_sm(freq))
#pylab.plot(freq,Eff_100_sm(freq))
#pylab.plot(freq,Eff_50_sm(freq)-dP*Eff_100_sm(freq))
#pylab.savefig(outdir+'June_'+date_ind+'_efftest',dpi=300)
#pylab.clf()


pylab.plot(freq,KA)
pylab.savefig(outdir+'June_'+date_ind+'_inttest',dpi=300)
pylab.clf()

pylab.plot(freq,Kt,'b--')
pylab.plot(freq,Kt50,'g')
pylab.plot(freq,Kt100,'r--')
pylab.xlim(50,100)
pylab.xlabel('Frequency (MHz)')
pylab.grid()
pylab.savefig(outdir+'June_'+date_ind+ '_caltest',dpi=300)
pylab.clf()

#pylab.plot(freq,(load_data - Eff_50_sm(freq)*KA)*Kt50)
#pylab.plot(freq,(term_data - 4*Eff_100_sm(freq)*KA)*Kt50)

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(load_full*Kt,mask=Kt_mask)),'b')
#ylab.plot(freq,Eff_50_sm(freq)*KA)
#ylab.plot(freq,Eff_100_sm(freq)*4*KA)
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(term_full*Kt,mask=Kt_mask)),'r')
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)), ma.compressed(ma.array((load_full - Eff_50_sm(freq)*KA)*Kt50,mask=Kt_mask)),'c')
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array((term_full-Eff_100_sm(freq)*4*KA)*Kt50,mask=Kt_mask)),'m--')
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array((short_full*Kt),mask=Kt_mask)),'g')
pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array((short_full*Kt50),mask=Kt_mask)),'y')
pylab.legend(('50 Omega','100 Omega','50 Omega w/ PZ','100 Omega w/ PZ','Short','Short w/ PZ'))
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(0,1100)
pylab.xlim(60,100)
pylab.xlabel('Frequency (MHz)')
pylab.grid()
pylab.savefig(outdir+'June_'+date_ind+'_comp_test',dpi=300)
pylab.clf()

fmin = where(freq<=60.)[0][-1]
fmax = where(freq<=100.)[0][-1]
old_avg = ma.mean(Kt[fmin:fmax]*short_full[fmin:fmax])
old_std = ma.std(Kt[fmin:fmax]*short_full[fmin:fmax])
new_avg = ma.mean(Kt50[fmin:fmax]*short_full[fmin:fmax])
new_std = ma.std(Kt50[fmin:fmax]*short_full[fmin:fmax])

mean_diff = ma.mean(Kt[fmin:fmax]*short_data[fmin:fmax]-Kt50[fmin:fmax]*short_data[fmin:fmax])
std_diff = ma.std(Kt[fmin:fmax]*short_data[fmin:fmax]-Kt50[fmin:fmax]*short_data[fmin:fmax])

print 'Old Avg Short is: ', old_avg, ' Kelvin'
print 'Old Short STD is: ', old_std, ' Kelvin'
print 'New Avg Short is: ', new_avg, ' Kelvin'
print 'New Short STD is: ', new_std, ' Kelvin'
print 'Avg Short Diff is: ', mean_diff, ' Kelvin'
print 'STD Short Diff is: ', std_diff, ' Kelvin'

f70 = where(freq<=70.)[0][-1]
print 'Old 70 MHz Short is: ', Kt[f70]*short_data[f70],' Kelvin'
print 'New 70 MHz Short is: ', Kt50[f70]*short_data[f70],' Kelvin'


#Kt_mask = ff.spike_flag(mean_Kt,mean_sm,freq,15.)
#print 'Percent of Kt mean data masked: ',100.*sum(Kt_mask)/len(Kt_mask)

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_sd*1.e9,mask=Kt_mask)))
pylab.xlim(60,120)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Power (nW)')
pylab.ylim(0,40)
pylab.grid()
#pylab.title('Time Mean without calibration')
pylab.savefig(outdir+'June_'+date_ind+'_mean_uncal_spectrum',dpi=300)
pylab.clf()

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt,mask=Kt_mask)))
pylab.xlim(60,120)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(0,6e3)
pylab.grid()
#pylab.title('Time Mean with JNC calibration')
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_spectrum',dpi=300)
pylab.clf()

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt,mask=Kt_mask)))
pylab.xlim(80,100)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(0.5e3,3e3)
pylab.grid()
#pylab.title('Time Mean with JNC calibration')
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_spectrum_zoom',dpi=300)
pylab.clf()

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_Kt,mask=Kt_mask)),c='b')
pylab.plot(freq,load_data*Kt,'g--')
pylab.plot(freq,short_data*Kt,'r--')
pylab.legend(('Sky','50 Ohm Load','Short'))
pylab.xlim(60,100)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(0,6e3)
pylab.grid()
#pylab.title('Time Mean with JNC calibration')
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_spectrum_ref',dpi=300)

pylab.ylim(0,1.1e4)
pylab.plot(freq,term_data*Kt,'c--')
pylab.plot(freq,noise_data*Kt,'m--')
pylab.legend(('Sky','50 Ohm Load','Short','100 Ohm Load','Noise Source')) 
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_spectrum_full_ref',dpi=300)
pylab.clf()

pylab.plot(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(mean_sd*1.e9,mask=Kt_mask)),c='b')
pylab.plot(freq,load_data*1.e9,'g--') 
pylab.plot(freq,short_data*1.e9,'r--') 
pylab.legend(('Sky','50 Ohm Load','Short')) 
pylab.xlim(60,100) 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Power (nW)') 
pylab.ylim(0,60) 
pylab.grid() 
#pylab.title('Time Mean with JNC calibration') 
pylab.savefig(outdir+'June_'+date_ind+'_mean_uncal_spectrum_ref',dpi=300)

pylab.ylim(0,60)
pylab.plot(freq,term_data*1.e9,'c--')
pylab.plot(freq,noise_data*1.e9,'m--')
pylab.legend(('Sky','50 Ohm Load','Short','100 Ohm Load','Noise Source'))
pylab.savefig(outdir+'June_'+date_ind+'_mean_uncal_spectrum_full_ref',dpi=300)

pylab.clf() 

pylab.plot(freq,load_data*1.e9,'g--') 
#pylab.scatter(freq,load_full*1.e9,c='g',edgecolor='g',s=3)
pylab.plot(freq,short_data*1.e9,'r--')
pylab.xlim(60,100) 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Power Fit (nW)') 
pylab.ylim(0,5) 
pylab.grid() 
pylab.plot(freq,term_data*1.e9,'c--')
pylab.legend(('50 Ohm Load','Short','100 Ohm Load'))
pylab.savefig(outdir+'June_'+date_ind+'_mean_uncal_ref_spectrum',dpi=300)

pylab.scatter(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(load_full,mask=Kt_mask))*1.e9,c='g',edgecolor='g',s=3) 
pylab.scatter(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(short_full,mask=Kt_mask))*1.e9,c='r',edgecolor='r',s=3) 
pylab.scatter(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(term_full,mask=Kt_mask))*1.e9,c='c',edgecolor='c',s=3) 
pylab.legend()
pylab.savefig(outdir+'June_'+date_ind+'_mean_uncal_ref_spectrum_data',dpi=300) 
pylab.clf()

pylab.plot(freq,load_data*Kt,'g--')  
pylab.plot(freq,short_data*Kt,'r--') 
pylab.xlim(60,100)  
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Temperature (Kelvin)')  
pylab.ylim(0,1000)   
pylab.grid()  
pylab.plot(freq,term_data*Kt,'c--')
pylab.legend(('50 Ohm Load','Short','100 Ohm Load'))
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_ref_spectrum',dpi=300)

pylab.scatter(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(load_full*Kt,mask=Kt_mask)),c='g',edgecolor='g',s=3)
pylab.scatter(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(short_full*Kt,mask=Kt_mask)),c='r',edgecolor='r',s=3)
pylab.scatter(ma.compressed(ma.array(freq,mask=Kt_mask)),ma.compressed(ma.array(term_full*Kt,mask=Kt_mask)),c='c',edgecolor='c',s=3)
pylab.legend() 
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_ref_spectrum_data',dpi=300)

pylab.clf()

pylab.plot(freq,load_data*Kt,'g--')
pylab.plot(freq,load_data*Kt50,'k')
pylab.plot(freq,short_data*Kt,'r--')
pylab.plot(freq,short_data*Kt50,'m')
pylab.xlim(60,100)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.ylim(0,1000)
pylab.grid()
pylab.plot(freq,term_data*Kt,'c--')
pylab.plot(freq,term_data*Kt50,'b')
#pylab.legend(('50 Ohm Load','Short','100 Ohm Load'))
pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_ref_spectrum_comp',dpi=300)
 
#pylab.scatter(freq,load_full*Kt,c='g',edgecolor='g',s=3)
#pylab.scatter(freq,short_full*Kt,c='r',edgecolor='r',s=3)
#pylab.scatter(freq,term_full*Kt,c='c',edgecolor='c',s=3)
#pylab.legend()
#pylab.savefig(outdir+'June_'+date_ind+'_mean_JNCcal_ref_spectrum_data',dpi=300)
 
pylab.clf()


fmin = where(freq<=60.)[0][-1]
fmax = where(freq<=120.)[0][-1]

pylab.imshow(sortdata[:,fmin:fmax]*1e9,vmin=0,vmax=60,aspect=60./(sorttime[-1]-sorttime[0]),extent=(60,120,sorttime[-1],sorttime[0]))
pylab.colorbar() 
pylab.title('Measured Antenna Power (nW) ') 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Sidereal Time (Hours)') 
pylab.savefig(outdir+'June_'+date_ind+'_unmasked_uncal_waterfall',dpi=300)
pylab.clf() 

 
fmin = where(freq<=70.)[0][-1]
pylab.scatter(ma.compressed(ma.array(sorttime,mask=sortmask[:,fmin])),ma.compressed(ma.array(sortdata[:,fmin]*1e9,mask=sortmask[:,fmin])),c='b',edgecolor='b',s=3)
pylab.xlim(0,24)
pylab.xlabel('Sidereal Time (Hours)')
pylab.ylim(0,60)
pylab.ylabel('Signal (nW)')
pylab.savefig(outdir+'June_'+date_ind+'_time_series_uncal_70mhz',dpi=300)
pylab.clf()


fmin = where(freq<=60.)[0][-1]
fmax = where(freq<=120.)[0][-1]
for i in range(0,len(sortind)):
    sortdata[i] = (sortdata[i]-short_data)*Kt*Eff_sm(freq)

pylab.imshow(sortdata[:,fmin:fmax],vmin=0,vmax=6e3,aspect=60./(sorttime[-1]-sorttime[0]),extent=(60,120,sorttime[-1],sorttime[0]))
pylab.colorbar()
pylab.title('Measured Sky Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.savefig(outdir+'June_'+date_ind+'_unmasked_cal_waterfall',dpi=300)
pylab.clf()

fmin2 = where(freq<=85.)[0][-1]
fmax2 = where(freq<=100.)[0][-1]
pylab.imshow(sortdata[:,fmin2:fmax2],vmin=1000,vmax=3e3,aspect=15./(sorttime[-1]-sorttime[0]),extent=(85,100,sorttime[-1],sorttime[0]))
pylab.colorbar()
pylab.title('Measured Sky Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.savefig(outdir+'June_'+date_ind+'_unmasked_cal_waterfall_FM',dpi=300)
pylab.clf()


for i in range(0,len(sortind)):
    for j in range(0,len(freq)):
        if sortmask[i,j]==1:
            sortdata[i,j] = 0

#pylab.imshow(sortdata[:,fmin:fmax]*1.e9,vmin=0,vmax=60,aspect=60./(sorttime[-1]-sorttime[0]),extent=(60,120,sorttime[-1],sorttime[0]))
#pylab.colorbar()
#pylab.title('Measured Antenna Power (nW), masked==0')
#pylab.xlabel('Frequency (MHz)')
#pylab.ylabel('Sidereal Time (Hours)')
#pylab.savefig(outdir+'June_'+date_ind+'_masked_uncal_waterfall',dpi=300)
#pylab.clf()

pylab.imshow(sortdata[:,fmin:fmax],vmin=0,vmax=6e3,aspect=60./(sorttime[-1]-sorttime[0]),extent=(60,120,sorttime[-1],sorttime[0]))
pylab.colorbar() 
pylab.title('Measured Sky Temperature (Kelvin), masked==0')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.savefig(outdir+'June_'+date_ind+'_masked_cal_waterfall',dpi=300)
pylab.clf()

fmin2 = where(freq<=85.)[0][-1]
fmax2 = where(freq<=100.)[0][-1]
pylab.imshow(sortdata[:,fmin2:fmax2],vmin=1000,vmax=3e3,aspect=15./(sorttime[-1]-sorttime[0]),extent=(85,100,sorttime[-1],sorttime[0]))
pylab.colorbar()
pylab.title('Measured Sky Temperature (Kelvin), masked==0')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.savefig(outdir+'June_'+date_ind+'_masked_waterfall_FM',dpi=300)
pylab.clf()

fmin3 = where(freq<=95.)[0][-1] 
fmax3 = where(freq<=115.)[0][-1] 
pylab.imshow(sortdata[:,fmin3:fmax3],vmin=1500,vmax=4.5e3,aspect=20./(sorttime[-1]-sorttime[0]),extent=(95,115,sorttime[-1],sorttime[0]))
pylab.colorbar() 
pylab.title('Measured Sky Temperature(Kelvin), masked==0')
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Sidereal Time (Hours)') 
pylab.savefig(outdir+'June_'+date_ind+'_masked_waterfall_RFI',dpi=300)
pylab.clf()

fmin = where(freq<=70.)[0][-1]
pylab.scatter(ma.compressed(ma.array(sorttime,mask=sortmask[:,fmin])),ma.compressed(ma.array(sortdata[:,fmin],mask=sortmask[:,fmin])),c='b',edgecolor='b',s=3)
pylab.xlim(0,24) 
pylab.xlabel('Sidereal Time (Hours)') 
pylab.ylim(0,6e3)
pylab.ylabel('Temperature (Kelvin)')
pylab.savefig(outdir+'June_'+date_ind+'_time_series_cal_70mhz',dpi=300)
pylab.clf() 

