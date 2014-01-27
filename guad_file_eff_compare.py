"""
Module to rebin a single day of data
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
#from pylab import *
#import matplotlib.pyplot as plt
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import sys
import ephem as eph
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly

caldir = '/lustre/tcv/truncated_data/'
maindir = '/lustre/tcv/freq_rebinned_data/'
newdir = '/lustre/tcv/mean_cal_data/'
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

timescale = 32
binscale = 32
directories = os.listdir(maindir)
date_ind = sys.argv[1]

#For efficiency calculation
#Try letting the efficiency phase delay float
#tdelay = arange(0,2e-9,1e-10)
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)

#GSM data (not needed here)
#gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_update.dat')
#gsm_freq = arange(60,90,1)
#gsm_time = arange(0,24,24./289)
#gsm_data = zeros((len(gsm_time),len(gsm_freq)))
#for t in range(0,len(gsm_time)):
#    for f in range(0,len(gsm_freq)):
#        gsm_data[t,f] = gsm_raw_data[f*289+t,2]

#Make data array
processed_data = []
processed_mask = []
processed_time = []
processed_volt = []
processed_data_b = []
processed_mask_b = []
processed_volt_b = []
processed_time_b = []
date_mid = [0.0,0.0,99.055,123.626,0.0,171.583,185.119,0.0,0.0,0.0,291.619,315.081,331.218,340.552]
split_short = [False,False,False,True,False,False,False,False,False,False,False,True,False,False]

if int(date_ind)<15:
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being analyzed is:',direct
    if split_short[int(date_ind)-1]:
        old_short_data = loadtxt(caldir+direct+'/mean_short_b.txt')
        old_short_data_b = loadtxt(caldir+direct+'/mean_short_a.txt')
        old_50ohm_data = loadtxt(caldir+direct+'/mean_50ohm_b.txt')
        old_50ohm_data_b = loadtxt(caldir+direct+'/mean_50ohm_a.txt')
    else:
        old_short_data = loadtxt(caldir+direct+'/mean_short.txt')
        old_50ohm_data = loadtxt(caldir+direct+'/mean_50ohm.txt')
    old_freq = arange(40,130,90./len(old_short_data))
    short_mask=zeros(len(old_short_data))
    short,smask,sfreq = fc.rebin(old_short_data,short_mask,old_freq,binscale)
    load,lmask,lfreq = fc.rebin(old_50ohm_data,short_mask,old_freq,binscale)
    if split_short[int(date_ind)-1]:    
        shortb,smaskb,sfreqb = fc.rebin(old_short_data_b,short_mask,old_freq,binscale)
        loadb,lmaskb,lfreqb = fc.rebin(old_50ohm_data_b,short_mask,old_freq,binscale)
    directory = maindir+direct+'/'
    dirlist = os.listdir(directory)
    for fname in dirlist:
        if fname.split('_')[-3] =='antenna':
            filename = directory+fname
            time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(filename)
            freq = arange(40.122075279756,129.9,0.244150559512)
            mask = fc.flagging(sub_data,freq,3.,binscale)  
            spike_mask = fc.spike_flag(sub_data,100.)
            for i in range(0,len(freq)):
                if spike_mask[i]==1.0:
                    mask[i] = 1.0
            if volt>10.0:
                if time>date_mid[int(date_ind)-1]:
                    processed_data.append(sub_data)
                    processed_mask.append(mask)
                    processed_time.append(time)
                    processed_volt.append(volt) 
                elif time<=date_mid[int(date_ind)-1]:
                    processed_data_b.append(sub_data)
                    processed_mask_b.append(mask)
                    processed_time_b.append(time)
                    processed_volt_b.append(volt)

print 'Number of datasets used is:',len(processed_time)

Kt = 300./(load-short)
if split_short[int(date_ind)-1]:
    Ktb = 300./(loadb-shortb)
#Kgu = loadtxt(maindir+direct+'/sm_un_sub_K.txt')
#Kgm = loadtxt(maindir+direct+'/sm_mean_sub_K.txt')

#Sidereal Time Calculation
sidereal_hour = []
idate = '2013/6/1'
for i in range(0,len(processed_time)):
    single_sid = fc.sidereal(processed_time[i],idate)
    sidereal_hour.append(single_sid)

sidereal_hour_b = []
for i in range(0,len(processed_time_b)):
    single_sid = fc.sidereal(processed_time_b[i],idate)
    sidereal_hour_b.append(single_sid)

#Sort data by LST
processed_data = array(processed_data)
processed_mask = array(processed_mask)
processed_volt = array(processed_volt)
processed_data_b = array(processed_data_b)
processed_mask_b = array(processed_mask_b)
processed_volt_b = array(processed_volt_b) 

sortind = argsort(sidereal_hour)
sorttime = zeros(len(processed_data))
sortdata = zeros((len(processed_data),len(processed_data[0])))
sortmask = zeros((len(processed_data),len(processed_data[0])))
sortvolt = zeros(len(processed_data))
sortindb = argsort(sidereal_hour_b)
sorttimeb = zeros(len(processed_data_b)) 
sortdatab = zeros((len(processed_data_b),len(processed_data[0])))
sortmaskb = zeros((len(processed_data_b),len(processed_data[0])))
sortvoltb = zeros(len(processed_data_b))
for i in range(0,len(sidereal_hour)):
    sorttime[i] = sidereal_hour[sortind[i]]
    sortdata[i] = processed_data[sortind[i]]
    sortmask[i] = processed_mask[sortind[i]]
    sortvolt[i] = processed_volt[sortind[i]]
for i in range(0,len(sidereal_hour_b)):
    sorttimeb[i] = sidereal_hour_b[sortindb[i]]
    sortdatab[i] = processed_data_b[sortindb[i]]
    sortmaskb[i] = processed_mask_b[sortindb[i]]
    sortvoltb[i] = processed_volt_b[sortindb[i]]

for i in range(1,len(sortvolt)):
    if (sortvolt[i]-sortvolt[i-1])>1.0:
        print 'Time index is ',i,' out of ',len(sortvolt)
        print 'Sidereal time of Battery Swap:', sorttime[i]
        print 'Raw time of Battery Swap:', processed_time[sortind[i]]

#Frequency Masking amount Calculation
percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency Masking',percent_masked

percent_maskedb = 100.*sum(sortmaskb)/(len(sortmaskb)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency Masking part2',percent_maskedb

#Adding in Time Masking
for i in range(0,len(freq)):
    new_mask = fc.timeflag(sortdata[:,i],sortmask[:,i],sorttime,3.,timescale)
    sortmask[:,i] = new_mask
    if date_mid[int(date_ind)-1]>0.0:
        new_maskb = fc.timeflag(sortdatab[:,i],sortmaskb[:,i],sorttimeb,3.,timescale) 
        sortmaskb[:,i] = new_maskb

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking',percent_masked 

percent_maskedb = 100.*sum(sortmaskb)/(len(sortmaskb)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking part2',percent_maskedb

#Calculating the time mean data:
mean_data =[]
for i in range(0,len(freq)):
    if sum(sortmask[:,i])==len(sortmask[:,i]):
        mean_data.append(0.0)
    else:
        single_data = ma.array(sortdata[:,i],mask=sortmask[:,i])
        single_comp = ma.compressed(single_data)
        single_mean = ma.mean(single_comp)
        mean_data.append(single_mean)

mean_datab =[]
for i in range(0,len(freq)): 
    if sum(sortmaskb[:,i])==len(sortmaskb[:,i]): 
        mean_datab.append(0.0) 
    else: 
        single_data = ma.array(sortdatab[:,i],mask=sortmaskb[:,i])
        single_comp = ma.compressed(single_data) 
        single_mean = ma.mean(single_comp) 
        mean_datab.append(single_mean)

#Masking the time mean data:
mean_data = array(mean_data)
mean_mask = zeros(len(mean_data))
spike_flag = fc.spike_flag(mean_data,5.)
for i in range(0,len(mean_data)):
    if mean_data[i] == 0.0:
        mean_mask[i] = 1.0
    if spike_flag[i] == 1.0:
        mean_mask[i] = 1.0
        print 'Spike Frequency at:',freq[i]
        mean_data[i] = 0.1*(mean_data[i+5]+mean_data[i+4]+mean_data[i+3]+mean_data[i+2]+mean_data[i+1]+mean_data[i-5]+mean_data[i-4]+mean_data[i-3]+mean_data[i-2]+mean_data[i-1])
print 'Number of Spiked Frequencies:',sum(spike_flag)
print 'Total Number of Masked Frequencies:',sum(mean_mask)

mean_datab = array(mean_datab) 
mean_maskb = zeros(len(mean_datab)) 
spike_flagb = fc.spike_flag(mean_datab,5.)
for i in range(0,len(mean_datab)):
    if mean_datab[i] == 0.0:
        mean_maskb[i] = 1.0
    if spike_flagb[i] == 1.0:
        mean_maskb[i] = 1.0
        print 'Spike Frequency in part 2 at:',freq[i]
        mean_datab[i] = 0.1*(mean_datab[i+5]+mean_datab[i+4]+mean_datab[i+3]+mean_datab[i+2]+mean_datab[i+1]+mean_datab[i-5]+mean_datab[i-4]+mean_datab[i-3]+mean_datab[i-2]+mean_datab[i-1])
print 'Number of Spiked Frequencies in part 2:',sum(spike_flagb) 
print 'Total Number of Masked Frequencies in part 2:',sum(mean_maskb)

def efficiency_test(td,ant_s11_file,amp_s_file):
    R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,td)
    R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
    Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
    Eff_sm = fc.smooth(Effic,F_ant,0.01)
    return Eff_sm

def residuals(td,ant_s11_file,amp_s_file,Kt,sh,freq,data):
    Eff_sm = efficiency_test(td,ant_s11_file,amp_s_file)
    mean_d = Kt*(data/Eff_sm(freq)-sh)
    minfreq = where(freq<=60)[0][-1]
    maxfreq = where(freq<=90)[0][-1]
    fit_d = poly.polyfit(log10(freq[minfreq:maxfreq]/70.),log10(mean_d[minfreq:maxfreq]),1)
    fit_val = 10**(poly.polyval(log10(freq/70.),fit_d))
    return mean_d[minfreq:maxfreq] - fit_val[minfreq:maxfreq]

bin = 9
mean_rebin,mask_rebin,new_freq = fc.rebin(mean_data,mean_mask,freq,bin)
mean_rebinb,mask_rebinb,new_freq = fc.rebin(mean_datab,mean_maskb,freq,bin)
short_rebin,smask_rebin,new_freq = fc.rebin(short,smask,freq,bin)
Kt_rebin,smask_rebin,new_freq = fc.rebin(Kt,smask,freq,bin)
if split_short[int(date_ind)-1]:
    short_rebinb,smask_rebinb,new_freq = fc.rebin(shortb,smaskb,freq,bin)
    Kt_rebinb,smask_rebinb,new_freq = fc.rebin(Ktb,smaskb,freq,bin)
freqs_rebin = new_freq
print freqs_rebin

td0 = [1e-12,]
test = opt.leastsq(residuals,td0[0],args=(ant_s11_file,amp_s_file,Kt_rebin,short_rebin,freqs_rebin,mean_rebin),maxfev=100000)
print test
if date_mid[int(date_ind)-1]>1.:
    if split_short[int(date_ind)-1]:
        testb = opt.leastsq(residuals,td0[0],args=(ant_s11_file,amp_s_file,Kt_rebinb,short_rebinb,freqs_rebin,mean_rebinb),maxfev=100000) 
    else:
        testb = opt.leastsq(residuals,td0[0],args=(ant_s11_file,amp_s_file,Kt_rebin,short_rebin,freqs_rebin,mean_rebinb),maxfev=100000)

    print testb


#print test



