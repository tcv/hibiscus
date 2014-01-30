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
#R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,3.129e-9)
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)
R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)

#GSM data (not needed here)
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_update.dat')
gsm_freq = arange(60,90,1)
gsm_time = arange(0,24,24./289)
gsm_data = zeros((len(gsm_time),len(gsm_freq)))
for t in range(0,len(gsm_time)):
    for f in range(0,len(gsm_freq)):
        gsm_data[t,f] = gsm_raw_data[f*289+t,2]

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
minfgu = [60.,70.,60.,60.,60.,60.,62.,60.,70.,62.,62.,60.,60.,60.]
maxfgu = [90.,90.,90.,90.,90.,90.,77.,90.,75.,90.,80.,90.,80.,90.]
minfgm = [60.,70.,60.,60.,65.,62.,62.,62.,70.,70.,70.,62.,62.,60.]
maxfgm = [90.,75.,75.,75.,70.,77.,77.,77.,75.,70.,75.,85.,77.,75.]
minft = [60.,70.,60.,60.,65.,62.,62.,62.,70.,62.,62.,62.,60.,62.]
maxft = [90.,90.,80.,80.,80.,77.,77.,80.,75.,80.,90.,82.,77.,80.]
minfgub = [70.,70.,60.,60.,70.,70.,60.,70.,70.,70.,70.,60.,60.,60.]
maxfgub = [70.,70.,80.,80.,70.,70.,90.,70.,70.,70.,82.,80.,80.,77.]
minfgmb = [70.,70.,60.,60.,70.,70.,62.,70.,70.,70.,70.,62.,62.,60.]
maxfgmb = [70.,70.,75.,75.,70.,70.,75.,70.,70.,70.,75.,80.,77.,77.]
minftb = [70.,70.,60.,60.,70.,70.,62.,70.,70.,70.,70.,62.,60.,62.] 
maxftb = [70.,70.,77.,75.,70.,70.,75.,70.,70.,70.,75.,82.,80.,75.] 

if int(date_ind)<15:
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being analyzed is:',direct
    if split_short[int(date_ind)-1]:
        old_short_data = loadtxt(caldir+direct+'/mean_short_b.txt')
        old_short_data_b = loadtxt(caldir+direct+'/mean_short_a.txt')
        old_50ohm_data = loadtxt(caldir+direct+'/mean_50ohm_b.txt')
        old_50ohm_data_b = loadtxt(caldir+direct+'/mean_50ohm_a.txt')
#        old_freq = arange(40,130,90./len(old_short_data_a))
    else:
        old_short_data = loadtxt(caldir+direct+'/mean_short.txt')
        old_50ohm_data = loadtxt(caldir+direct+'/mean_50ohm.txt')
#    old_short_data_a = loadtxt(caldir+direct+'/mean_short_a.txt')
#    old_short_data_b = loadtxt(caldir+direct+'/mean_short_b.txt')
    old_freq = arange(40,130,90./len(old_short_data))
    short_mask=zeros(len(old_short_data))
    short,smask,sfreq = fc.rebin(old_short_data,short_mask,old_freq,binscale)
    load,lmask,lfreq = fc.rebin(old_50ohm_data,short_mask,old_freq,binscale)
    if split_short[int(date_ind)-1]:    
        shortb,smaskb,sfreqb = fc.rebin(old_short_data_b,short_mask,old_freq,binscale)
        loadb,lmaskb,lfreqb = fc.rebin(old_50ohm_data_b,short_mask,old_freq,binscale)
    directory = maindir+direct+'/'
#    new_directory = newdir+direct+'/'
    dirlist = os.listdir(directory)
    for fname in dirlist:
        if fname.split('_')[-3] =='antenna':
            filename = directory+fname
            time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(filename)
#            width = 90.0/len(sub_data)
#            freq = arange(40,130.0,width)
            freq = arange(40.122075279756,129.9,0.244150559512)
            mask = fc.flagging(sub_data,freq,3.,binscale)  
            if split_short[int(date_ind)-1]:
                if time>date_mid[int(date_ind)-1]:
                    new_data = (sub_data/Eff_sm(freq)-shortb)
                else: 
                    new_data = (sub_data/Eff_sm(freq)-short)
            else:
                new_data = (sub_data/Eff_sm(freq)-short)
            spike_mask = fc.spike_flag(new_data,100.)
            for i in range(0,len(freq)):
                if spike_mask[i]==1.0:
                    mask[i] = 1.0
            if volt>10.0:
                if time>date_mid[int(date_ind)-1]:
                    processed_data.append(new_data)
                    processed_mask.append(mask)
                    processed_time.append(time)
                    processed_volt.append(volt) 
                elif time<=date_mid[int(date_ind)-1]:
                    processed_data_b.append(new_data)
                    processed_mask_b.append(mask)
                    processed_time_b.append(time)
                    processed_volt_b.append(volt)

print 'Number of datasets used is:',len(processed_time)

Kt = 300./(load-short)
if split_short[int(date_ind)-1]:
    Ktb = 300./(loadb-shortb)
Kgu = loadtxt(maindir+direct+'/sm_un_sub_K.txt')
Kgm = loadtxt(maindir+direct+'/sm_mean_sub_K.txt')
#print Kgm
Kgm = abs(Kgm)

#print where(isnan(Kgm))
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
#if date_mid[int(date_ind)-1]>0.0:
for i in range(0,len(freq)):
    new_mask = fc.timeflag(sortdata[:,i],sortmask[:,i],sorttime,2.5,timescale)
    sortmask[:,i] = new_mask
    if date_mid[int(date_ind)-1]>0.0:
        new_maskb = fc.timeflag(sortdatab[:,i],sortmaskb[:,i],sorttimeb,2.5,timescale) 
        sortmaskb[:,i] = new_maskb

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking',percent_masked 

percent_maskedb = 100.*sum(sortmaskb)/(len(sortmaskb)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking part2',percent_maskedb

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

mean_data = array(mean_data)
mean_mask = zeros(len(mean_data))
spike_flag = fc.spike_flag(mean_data,4.)
for i in range(0,len(mean_data)):
    if mean_data[i] == 0.0:
        mean_mask[i] = 1.0
    if spike_flag[i] == 1.0:
        mean_mask[i] = 1.0
        print 'Spike Frequency at:',freq[i]
        mean_data[i] = 0.1*(mean_data[i+5]+mean_data[i+4]+mean_data[i+3]+mean_data[i+2]+mean_data[i+1]+mean_data[i-5]+mean_data[i-4]+mean_data[i-3]+mean_data[i-2]+mean_data[i-1])
#spike_flag = fc.spike_flag(mean_data,100)
print 'Number of Spiked Frequencies:',sum(spike_flag)
print 'Total Number of Masked Frequencies:',sum(mean_mask)

mean_datab = array(mean_datab) 
mean_maskb = zeros(len(mean_datab)) 
spike_flagb = fc.spike_flag(mean_datab,4.)
for i in range(0,len(mean_datab)):
    if mean_datab[i] == 0.0:
        mean_maskb[i] = 1.0
    if spike_flagb[i] == 1.0:
        mean_maskb[i] = 1.0
        print 'Spike Frequency in part 2 at:',freq[i]
        mean_datab[i] = 0.1*(mean_datab[i+5]+mean_datab[i+4]+mean_datab[i+3]+mean_datab[i+2]+mean_datab[i+1]+mean_datab[i-5]+mean_datab[i-4]+mean_datab[i-3]+mean_datab[i-2]+mean_datab[i-1])
#spike_flag = fc.spike_flag(mean_data,100)
print 'Number of Spiked Frequencies in part 2:',sum(spike_flagb) 
print 'Total Number of Masked Frequencies in part 2:',sum(mean_maskb)

if split_short[int(date_ind)-1]:
    tcal_mean = mean_data*Kt
    tcal_meanb = mean_datab*Ktb
else: 
    tcal_mean = mean_data*Kt
    tcal_meanb = mean_datab*Kt
gucal_mean = mean_data*Kgu
gmcal_mean = mean_data*Kgm
gucal_meanb = mean_datab*Kgu
gmcal_meanb = mean_datab*Kgm

bin = 8
tcal_mean_rebin,tcal_mask_rebin,new_freq = fc.rebin(tcal_mean,mean_mask,freq,bin)
gucal_mean_rebin,gucal_mask_rebin,new_freq = fc.rebin(gucal_mean,mean_mask,freq,bin)
gmcal_mean_rebin,gmcal_mask_rebin,new_freq = fc.rebin(gmcal_mean,mean_mask,freq,bin)
tcal_mean_rebinb,tcal_mask_rebinb,new_freq = fc.rebin(tcal_meanb,mean_maskb,freq,bin)
gucal_mean_rebinb,gucal_mask_rebinb,new_freq = fc.rebin(gucal_meanb,mean_maskb,freq,bin)
gmcal_mean_rebinb,gmcal_mask_rebinb,new_freq = fc.rebin(gmcal_meanb,mean_maskb,freq,bin)
freqs_rebin = new_freq
print freqs_rebin


savetxt(newdir+direct+'_tcal_mean.txt',tcal_mean_rebin,delimiter=' ')
savetxt(newdir+direct+'_gucal_mean.txt',gucal_mean_rebin,delimiter=' ')
savetxt(newdir+direct+'_gmcal_mean.txt',gmcal_mean_rebin,delimiter=' ')
if date_mid[int(date_ind)-1]>0.0:
    savetxt(newdir+direct+'_tcal_meanb.txt',tcal_mean_rebinb,delimiter=' ')
    savetxt(newdir+direct+'_gucal_meanb.txt',gucal_mean_rebinb,delimiter=' ')
    savetxt(newdir+direct+'_gmcal_meanb.txt',gmcal_mean_rebinb,delimiter=' ')

minfreq = where(freqs_rebin<=minft[int(date_ind)-1])[0][-1]
maxfreq = where(freqs_rebin<=maxft[int(date_ind)-1])[0][-1]
minfreqb = where(freqs_rebin<=minftb[int(date_ind)-1])[0][-1]
maxfreqb = where(freqs_rebin<=maxftb[int(date_ind)-1])[0][-1]

if maxfreq>minfreq:
    tcal_fit1 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(tcal_mean_rebin[minfreq:maxfreq]),1)
    tcal_fit2 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(tcal_mean_rebin[minfreq:maxfreq]),2)
    tcal_fit3 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(tcal_mean_rebin[minfreq:maxfreq]),3)
    tcal_fit4 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(tcal_mean_rebin[minfreq:maxfreq]),4)
    tcal_fit5 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(tcal_mean_rebin[minfreq:maxfreq]),5)
    tcal_fit6 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(tcal_mean_rebin[minfreq:maxfreq]),6)
    print 'Johnson Noise Cal Fit Params are:',tcal_fit1
    print 'Johnson Noise Cal Fit Params are:',tcal_fit2
    print 'Johnson Noise Cal Fit Params are:',tcal_fit3
    print 'Johnson Noise Cal Fit Params are:',tcal_fit4
    print 'Johnson Noise Cal Fit Params are:',tcal_fit5
    print 'Johnson Noise Cal Fit Params are:',tcal_fit6

if maxfreqb>minfreqb:
    if date_mid[int(date_ind)-1]>0.0:
        tcal_fit1b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(tcal_mean_rebinb[minfreqb:maxfreqb]),1)
        tcal_fit2b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(tcal_mean_rebinb[minfreqb:maxfreqb]),2)
        tcal_fit3b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(tcal_mean_rebinb[minfreqb:maxfreqb]),3)
        tcal_fit4b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(tcal_mean_rebinb[minfreqb:maxfreqb]),4)
        tcal_fit5b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(tcal_mean_rebinb[minfreqb:maxfreqb]),5)
        tcal_fit6b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(tcal_mean_rebinb[minfreqb:maxfreqb]),6)
        print 'Johnson Noise Cal Fit Params (part 2) are:',tcal_fit1b
        print 'Johnson Noise Cal Fit Params (part 2) are:',tcal_fit2b
        print 'Johnson Noise Cal Fit Params (part 2) are:',tcal_fit3b
        print 'Johnson Noise Cal Fit Params (part 2) are:',tcal_fit4b
        print 'Johnson Noise Cal Fit Params (part 2) are:',tcal_fit5b
        print 'Johnson Noise Cal Fit Params (part 2) are:',tcal_fit6b

if maxfreq==minfreq:
    tcal_fit1 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),1)
    tcal_fit2 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),2)
    tcal_fit3 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),3)
    tcal_fit4 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),4)
    tcal_fit5 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),5)
    tcal_fit6 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),6)
if maxfreqb==minfreqb:
    tcal_fit1b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),1) 
    tcal_fit2b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),2)
    tcal_fit3b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),3)
    tcal_fit4b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),4)
    tcal_fit5b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),5)
    tcal_fit6b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),6)

minfreq = where(freqs_rebin<=minfgu[int(date_ind)-1])[0][-1]
maxfreq = where(freqs_rebin<=maxfgu[int(date_ind)-1])[0][-1] 
minfreqb = where(freqs_rebin<=minfgub[int(date_ind)-1])[0][-1] 
maxfreqb = where(freqs_rebin<=maxfgub[int(date_ind)-1])[0][-1] 

if maxfreq>minfreq:
    gucal_fit1 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gucal_mean_rebin[minfreq:maxfreq]),1)
    gucal_fit2 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gucal_mean_rebin[minfreq:maxfreq]),2)
    gucal_fit3 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gucal_mean_rebin[minfreq:maxfreq]),3)
    gucal_fit4 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gucal_mean_rebin[minfreq:maxfreq]),4)
    gucal_fit5 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gucal_mean_rebin[minfreq:maxfreq]),5)
    gucal_fit6 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gucal_mean_rebin[minfreq:maxfreq]),6)
    print 'GSM without Mean Sub Cal Fit Params are:',gucal_fit1
    print 'GSM without Mean Sub Cal Fit Params are:',gucal_fit2
    print 'GSM without Mean Sub Cal Fit Params are:',gucal_fit3
    print 'GSM without Mean Sub Cal Fit Params are:',gucal_fit4
    print 'GSM without Mean Sub Cal Fit Params are:',gucal_fit5
    print 'GSM without Mean Sub Cal Fit Params are:',gucal_fit6

if maxfreqb>minfreqb:
    if date_mid[int(date_ind)-1]>0.0:
        gucal_fit1b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gucal_mean_rebinb[minfreqb:maxfreqb]),1)
        gucal_fit2b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gucal_mean_rebinb[minfreqb:maxfreqb]),2)
        gucal_fit3b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gucal_mean_rebinb[minfreqb:maxfreqb]),3)
        gucal_fit4b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gucal_mean_rebinb[minfreqb:maxfreqb]),4)
        gucal_fit5b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gucal_mean_rebinb[minfreqb:maxfreqb]),5)
        gucal_fit6b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gucal_mean_rebinb[minfreqb:maxfreqb]),6)
        print 'GSM without Mean Sub Cal Fit Params (part 2) are:',gucal_fit1b
        print 'GSM without Mean Sub Cal Fit Params (part 2) are:',gucal_fit2b
        print 'GSM without Mean Sub Cal Fit Params (part 2) are:',gucal_fit3b
        print 'GSM without Mean Sub Cal Fit Params (part 2) are:',gucal_fit4b
        print 'GSM without Mean Sub Cal Fit Params (part 2) are:',gucal_fit5b
        print 'GSM without Mean Sub Cal Fit Params (part 2) are:',gucal_fit6b

if maxfreq==minfreq:
    gucal_fit1 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),1) 
    gucal_fit2 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),2) 
    gucal_fit3 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),3) 
    gucal_fit4 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),4) 
    gucal_fit5 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),5) 
    gucal_fit6 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),6) 
if maxfreqb==minfreqb: 
    gucal_fit1b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),1)
    gucal_fit2b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),2) 
    gucal_fit3b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),3) 
    gucal_fit4b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),4) 
    gucal_fit5b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),5) 
    gucal_fit6b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),6) 

minfreq = where(freqs_rebin<=minfgm[int(date_ind)-1])[0][-1] 
maxfreq = where(freqs_rebin<=maxfgm[int(date_ind)-1])[0][-1]  
minfreqb = where(freqs_rebin<=minfgmb[int(date_ind)-1])[0][-1]  
maxfreqb = where(freqs_rebin<=maxfgmb[int(date_ind)-1])[0][-1]  

if maxfreq>minfreq:
    gmcal_fit1 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gmcal_mean_rebin[minfreq:maxfreq]),1)
    gmcal_fit2 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gmcal_mean_rebin[minfreq:maxfreq]),2)
    gmcal_fit3 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gmcal_mean_rebin[minfreq:maxfreq]),3)
    gmcal_fit4 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gmcal_mean_rebin[minfreq:maxfreq]),4)
    gmcal_fit5 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gmcal_mean_rebin[minfreq:maxfreq]),5)
    gmcal_fit6 = poly.polyfit(log10(freqs_rebin[minfreq:maxfreq]/70.),log10(gmcal_mean_rebin[minfreq:maxfreq]),6)
    print 'GSM with Mean Sub Cal Fit Params are:',gmcal_fit1
    print 'GSM with Mean Sub Cal Fit Params are:',gmcal_fit2
    print 'GSM with Mean Sub Cal Fit Params are:',gmcal_fit3
    print 'GSM with Mean Sub Cal Fit Params are:',gmcal_fit4
    print 'GSM with Mean Sub Cal Fit Params are:',gmcal_fit5
    print 'GSM with Mean Sub Cal Fit Params are:',gmcal_fit6

if maxfreqb>minfreqb:
    if date_mid[int(date_ind)-1]>0.0:
        gmcal_fit1b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gmcal_mean_rebinb[minfreqb:maxfreqb]),1)
        gmcal_fit2b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gmcal_mean_rebinb[minfreqb:maxfreqb]),2)
        gmcal_fit3b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gmcal_mean_rebinb[minfreqb:maxfreqb]),3)
        gmcal_fit4b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gmcal_mean_rebinb[minfreqb:maxfreqb]),4)
        gmcal_fit5b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gmcal_mean_rebinb[minfreqb:maxfreqb]),5)
        gmcal_fit6b = poly.polyfit(log10(freqs_rebin[minfreqb:maxfreqb]/70.),log10(gmcal_mean_rebinb[minfreqb:maxfreqb]),6)
        print 'GSM with Mean Sub Cal Fit Params (part 2) are:',gmcal_fit1b
        print 'GSM with Mean Sub Cal Fit Params (part 2) are:',gmcal_fit2b
        print 'GSM with Mean Sub Cal Fit Params (part 2) are:',gmcal_fit3b
        print 'GSM with Mean Sub Cal Fit Params (part 2) are:',gmcal_fit4b
        print 'GSM with Mean Sub Cal Fit Params (part 2) are:',gmcal_fit5b
        print 'GSM with Mean Sub Cal Fit Params (part 2) are:',gmcal_fit6b

if maxfreq==minfreq: 
    gmcal_fit1 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),1)
    gmcal_fit2 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),2)
    gmcal_fit3 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),3)
    gmcal_fit4 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),4)
    gmcal_fit5 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),5)
    gmcal_fit6 = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),6)
if maxfreqb==minfreqb:
    gmcal_fit1b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),1)
    gmcal_fit2b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),2)
    gmcal_fit3b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),3)
    gmcal_fit4b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),4)
    gmcal_fit5b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),5)
    gmcal_fit6b = poly.polyfit(log10(freqs_rebin/10.),log10(ones(len(freqs_rebin))),6)

savetxt(newdir+direct+'_tcal_res1.txt',tcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit1)),delimiter=' ')
savetxt(newdir+direct+'_tcal_res2.txt',tcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit2)),delimiter=' ')
savetxt(newdir+direct+'_tcal_res3.txt',tcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit3)),delimiter=' ')
savetxt(newdir+direct+'_tcal_res4.txt',tcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit4)),delimiter=' ')
savetxt(newdir+direct+'_tcal_res5.txt',tcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit5)),delimiter=' ')
savetxt(newdir+direct+'_tcal_res6.txt',tcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit6)),delimiter=' ')

savetxt(newdir+direct+'_gucal_res1.txt',gucal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit1)),delimiter=' ')
savetxt(newdir+direct+'_gucal_res2.txt',gucal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit2)),delimiter=' ')
savetxt(newdir+direct+'_gucal_res3.txt',gucal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit3)),delimiter=' ')
savetxt(newdir+direct+'_gucal_res4.txt',gucal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit4)),delimiter=' ')
savetxt(newdir+direct+'_gucal_res5.txt',gucal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit5)),delimiter=' ')
savetxt(newdir+direct+'_gucal_res6.txt',gucal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit6)),delimiter=' ')

savetxt(newdir+direct+'_gmcal_res1.txt',gmcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit1)),delimiter=' ')
savetxt(newdir+direct+'_gmcal_res2.txt',gmcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit2)),delimiter=' ')
savetxt(newdir+direct+'_gmcal_res3.txt',gmcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit3)),delimiter=' ')
savetxt(newdir+direct+'_gmcal_res4.txt',gmcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit4)),delimiter=' ')
savetxt(newdir+direct+'_gmcal_res5.txt',gmcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit5)),delimiter=' ')
savetxt(newdir+direct+'_gmcal_res6.txt',gmcal_mean_rebin-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit6)),delimiter=' ')

if date_mid[int(date_ind)-1]>0.0:
    savetxt(newdir+direct+'_tcal_res1b.txt',tcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit1b)),delimiter=' ')
    savetxt(newdir+direct+'_tcal_res2b.txt',tcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit2b)),delimiter=' ')
    savetxt(newdir+direct+'_tcal_res3b.txt',tcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit3b)),delimiter=' ')
    savetxt(newdir+direct+'_tcal_res4b.txt',tcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit4b)),delimiter=' ')
    savetxt(newdir+direct+'_tcal_res5b.txt',tcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit5b)),delimiter=' ')
    savetxt(newdir+direct+'_tcal_res6b.txt',tcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),tcal_fit6b)),delimiter=' ')

    savetxt(newdir+direct+'_gucal_res1b.txt',gucal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit1b)),delimiter=' ')
    savetxt(newdir+direct+'_gucal_res2b.txt',gucal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit2b)),delimiter=' ')
    savetxt(newdir+direct+'_gucal_res3b.txt',gucal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit3b)),delimiter=' ')
    savetxt(newdir+direct+'_gucal_res4b.txt',gucal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit4b)),delimiter=' ')
    savetxt(newdir+direct+'_gucal_res5b.txt',gucal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit5b)),delimiter=' ')
    savetxt(newdir+direct+'_gucal_res6b.txt',gucal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gucal_fit6b)),delimiter=' ')
 
    savetxt(newdir+direct+'_gmcal_res1b.txt',gmcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit1b)),delimiter=' ')
    savetxt(newdir+direct+'_gmcal_res2b.txt',gmcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit2b)),delimiter=' ')
    savetxt(newdir+direct+'_gmcal_res3b.txt',gmcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit3b)),delimiter=' ')
    savetxt(newdir+direct+'_gmcal_res4b.txt',gmcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit4b)),delimiter=' ')
    savetxt(newdir+direct+'_gmcal_res5b.txt',gmcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit5b)),delimiter=' ')
    savetxt(newdir+direct+'_gmcal_res6b.txt',gmcal_mean_rebinb-10**(poly.polyval(log10(freqs_rebin/70.),gmcal_fit6b)),delimiter=' ')

