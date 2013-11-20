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

caldir = '/lustre/tcv/truncated_data/'
maindir = '/lustre/tcv/freq_rebinned_data/'
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

timescale = 32
binscale = 32
directories = os.listdir(maindir)
date_ind = sys.argv[1]

# For efficiency calculation
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)
R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)

#GSM data rearranged to array with first index = time, second index = freq
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_update.dat')
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
if int(date_ind)<15:
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being analyzed is:',direct
    old_short_data = loadtxt(caldir+direct+'/mean_short.txt')
    old_freq = arange(40,130,90./len(old_short_data))
    short_mask=zeros(len(old_short_data))
    short,smask,sfreq = fc.rebin(old_short_data,short_mask,old_freq,binscale)
    directory = maindir+direct+'/'
    dirlist = os.listdir(directory)
    for fname in dirlist:
        if fname.split('_')[-3] =='antenna':
            filename = directory+fname
            time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(filename)
#            width = 90.0/len(sub_data)
#            freq = arange(40,130.0,width)
            freq = arange(40.122075279756,129.9,0.244150559512)
#            print shape(sub_data),shape(freq)
            mask = fc.flagging(sub_data,freq,3.,binscale)
            new_data = (sub_data/Eff_sm(freq)-short)
            spike_mask = fc.spike_flag(new_data,50.)
            for i in range(0,len(freq)):
                if spike_mask[i]==1.0:
                    mask[i] = 1.0
            if volt>10.0:
                processed_data.append(new_data)
                processed_mask.append(mask)
                processed_time.append(time)

#Sidereal Time Calculation
initial = eph.date('2013/6/1')
guad = eph.Observer()
guad.lon = '-118.3'
guad.lat = '28.8833'

sidereal_hour = []
for i in range(0,len(processed_time)):
    single_date = eph.date(initial+processed_time[i]/24.)
    guad.date = single_date
    single_time = guad.sidereal_time()
    sidereal_hour.append(single_time*12./pi)

#Sort data by LST
processed_data = array(processed_data)
processed_mask = array(processed_mask)

sortind = argsort(sidereal_hour)
sorttime = zeros(len(processed_data))
sortdata = zeros((len(processed_data),len(processed_data[0])))
sortmask = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(sidereal_hour)):
    sorttime[i] = sidereal_hour[sortind[i]]
    sortdata[i] = processed_data[sortind[i]]
    sortmask[i] = processed_mask[sortind[i]]

#print where(isnan(sortdata))
#Frequency Masking Amount Calculation
percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency Masking',percent_masked

#Adding in Time Masking
for i in range(0,len(freq)):
    new_mask = fc.timeflag(sortdata[:,i],sortmask[:,i],sorttime,3.,timescale)
    sortmask[:,i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking',percent_masked 

#Time rebinning to match gsm binning (5 minutes)
stack_time = gsm_time
stack_data = zeros((len(stack_time),len(freq)))
for i in range(0,len(stack_time)):
    sub_data = []
    sub_mask = []
    num_mean = 0
    for j in range(0,len(sorttime)):
       if abs(stack_time[i]-sorttime[j])<=(stack_time[1]-stack_time[0])/2.:
           sub_data.append(sortdata[j])
           sub_mask.append(sortmask[j])
           num_mean = num_mean+1
    sub_data = array(sub_data)
    sub_mask = array(sub_mask)
    if num_mean>=1.:
        for f in range(0,len(freq)):
            if sum(sub_mask[:,f])==len(sub_mask[:,f]):
                stack_data[i,f] = 0.0
            else: 
                single_data = ma.array(sub_data[:,f],mask=sub_mask[:,f])
                single_comp = ma.compressed(single_data)
                stack_data[i,f] = ma.mean(single_comp)

#Only grabbing the gsm times that we have measured data for:
lim_stack = []
lim_gsm = []
#lim_stack_mask = []
for i in range(0,len(stack_data)):
    if sum(stack_data[i])>0: 
        lim_stack.append(stack_data[i])
#        single_mask = fc.spike_flag(stack_data[i],50.)
#        lim_stack_mask.append(single_mask)
        single_smooth = itp.UnivariateSpline(gsm_freq,gsm_data[i])
        lim_gsm.append(single_smooth(freq))
lim_stack = array(lim_stack) 
lim_gsm = array(lim_gsm)

stack_mask = zeros((len(lim_stack),len(lim_stack[0])))
for i in range(0,len(stack_mask)):
    for j in range(0,len(stack_mask[0])):
        if lim_stack[i,j]<=0.0:
            stack_mask[i,j] = 1.0

#print sum(stack_mask)
mean_stack_data = []
mean_gsm = []
for i in range(0,len(lim_stack[0])):
    if len(stack_mask[:,i])==sum(stack_mask[:,i]):
        mean_stack_data.append(0.0)
        mean_gsm.append(0.0)
    else:
        single_lim_stack = ma.array(lim_stack[:,i],mask=stack_mask[:,i])
        single_lim_comp = ma.compressed(single_lim_stack)
        mean_stack_data.append(ma.mean(single_lim_comp))
        single_lim_gsm = ma.array(lim_gsm[:,i],mask=stack_mask[:,i])
        single_gsm_comp = ma.compressed(single_lim_gsm)
        mean_gsm.append(ma.mean(single_gsm_comp))
mean_stack_data = array(mean_stack_data)
mean_gsm = array(mean_gsm)
#mean_stack_data = ma.mean(lim_stack,axis=0)

pylab.plot(freq,mean_stack_data*10**9)
pylab.grid()
pylab.xlim(40,130)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Raw Data (nW)')
pylab.title('Time Averaged Power Spectrum')
pylab.savefig(maindir+direct+'_mean_spectra',dpi=300)
pylab.clf()

pylab.plot(freq,mean_gsm)
pylab.grid()
pylab.xlim(40,130)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (K)')
pylab.title('Time Averaged GSM Signal')
pylab.savefig(maindir+direct+'_mean_gsm_spectra',dpi=300)
pylab.clf()

#print shape(stack_data),shape(mean_stack_data)
lim_ms_stack= zeros((len(lim_stack),len(stack_data[0])))
for i in range(0,len(lim_stack)):
    lim_ms_stack[i] = lim_stack[i]-mean_stack_data

#mean_gsm = ma.mean(lim_gsm,axis=0)
#mean_gsm_smooth = itp.UnivariateSpline(gsm_freq,mean_gsm)
lim_ms_gsm = zeros((len(lim_stack),len(freq)))
for i in range(0,len(lim_stack)):
    lim_ms_gsm[i] = lim_gsm[i]-mean_gsm

lim_ms_stack = array(lim_ms_stack)
lim_ms_gsm = array(lim_ms_gsm)

unsubK = []
meansubK = []
Ko0 = [1e10]
Km0 = [1e10]
for f in range(0,len(freq)):
    fKo = lambda Ko,d,g: Ko*d-g
    fKm = lambda Km,d,g: Km*d-g
#    Ko0 = [1e10]
#    Km0 = [1e10]
    if len(stack_mask[:,f])==sum(stack_mask[:,f]):
        unsubK.append(Ko0)
        meansubK.append(Km0)
        print freq[f]
    else:        
        do_array = ma.array(lim_stack[:,f],mask=stack_mask[:,f])
        do = ma.compressed(do_array)
        dm_array = ma.array(lim_ms_stack[:,f],mask=stack_mask[:,f])
        dm = ma.compressed(dm_array)
        go_array = ma.array(lim_gsm[:,f],mask=stack_mask[:,f])
        go = ma.compressed(go_array)
        gm_array = ma.array(lim_ms_gsm[:,f],mask=stack_mask[:,f])
        gm = ma.compressed(gm_array)
        Ko = opt.leastsq(fKo,Ko0,args=(do,go),maxfev=100000)
        Km = opt.leastsq(fKm,Km0,args=(dm,gm),maxfev=100000)
        Ko0 = Ko[0]
        Km0 = Km[0]
        unsubK.append(Ko[0])
        meansubK.append(Km[0])
    
unsubK = array(unsubK)
meansubK = array(meansubK)
f70 = where(freq<=70.)[0][-1]
print 'unsubK at 70 MHz is:',unsubK[f70]
print 'meansubK at 70 MHz is:',meansubK[f70]

savetxt(directory+'un_sub_K.txt',unsubK,delimiter = ' ')
savetxt(directory+'mean_sub_K.txt',meansubK,delimiter=' ')

