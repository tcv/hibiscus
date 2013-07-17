from numpy import *
import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf

maindir = 'Isla_Guadalupe_data_jun_2013/June03_day_to_night/'
maindir2 = 'Isla_Guadalupe_data_jun_2013/June04_day_to_night/'
rebin_ant = []
ant_time = []
ant_mask = []
load = []
load_time = []
term = []
term_time = []
short = []
short_time = []
noise = []
noise_time = []
binscale = 20
timescale = 10
ant_limit = []
mask_limit = []
time_limit = []
dataset = 0

directories = os.listdir(maindir)
for direct in directories:
    if direct.split('-')[0]=='2013':
        directory = maindir+direct+'/'
        print directory
        print shape(rebin_ant)
        if len(direct.split('_'))<2:
            dirlist = os.listdir(directory)
            for fname in dirlist:
                filename=directory+fname
                time,form,sub_data,mask,freq,volt = fc.loadsingle(filename)
                width = 250.0/len(sub_data)
                freqs = arange(0,250.0,width)
                if len(sub_data)>1:
                    if form=='antenna':
#                mask = zeros(len(sub_data))
                        masked_sub = 10**(sub_data/10.)
                        mask = fc.flagging(masked_sub,freqs,3.,30)
                        excess_mask = fc.spike_flag(new_data,100)
                        for i in range(0,len(excess_mask)):
                            if excess_mask[i]==1.0:
                                mask[i] = 1.0
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        if len(time_limit)<1:
                            time_test = time
                        else:
                            time_test = time_limit[-1]
                        if time-time_test>0.1:
                            timebin_data,timebin_mask = fc.timerebin(ant_limit,mask_limit)
                            rebin_ant.append(timebin_data)
                            ant_time.append(ma.mean(time_limit[0:-2]))
                            ant_mask.append(timebin_mask)
                            ant_limit = []
                            mask_limit = []
                            time_limit = []
                            dataset = 0                            
                        time_limit.append(time)
                        ant_limit.append(new_data)
                        mask_limit.append(new_mask)
                        dataset = dataset + 1
                        if dataset == timescale:
                            timebin_data,timebin_mask = fc.timerebin(ant_limit,mask_limit)
                            rebin_ant.append(timebin_data)
                            ant_time.append(ma.mean(time_limit))
                            ant_mask.append(timebin_mask)
                            ant_limit = []
                            mask_limit = []
                            time_limit = []
                            dataset = 0
                    elif form=='50ohm':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        load.append(new_data)
                        load_time.append(time)
                    elif form=='open':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        term.append(new_data)
                        term_time.append(time)
                    elif form=='short':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        short.append(new_data)
                        short_time.append(time)
                    elif form=='noise':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        noise.append(new_data)
                        noise_time.append(time)
                        
directories2 = os.listdir(maindir2)
for direct in directories2:
    if direct.split('-')[0]=='2013':
        directory = maindir2+direct+'/'
        print directory
        print shape(rebin_ant)
        if len(direct.split('_'))<2:
            dirlist = os.listdir(directory)
            for fname in dirlist:
                filename=directory+fname
                time,form,sub_data,mask,freq,volt = fc.loadsingle(filename)
                width = 250.0/len(sub_data)
                freqs = arange(0,250.0,width)
                if len(sub_data)>1:
                    if form=='antenna':
#                mask = zeros(len(sub_data))
                        masked_sub = 10**(sub_data/10.)
                        mask = fc.flagging(masked_sub,freqs,3.,30)
                        excess_mask = fc.spike_flag(new_data,25)
                        for i in range(0,len(excess_mask)):
                            if excess_mask[i]==1.0:
                                mask[i] = 1.0
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        if len(time_limit)<1:
                            time_test = time
                        else:
                            time_test = time_limit[-1]
                        if time-time_test>0.1:
                            timebin_data,timebin_mask = fc.timerebin(ant_limit,mask_limit)
                            rebin_ant.append(timebin_data)
                            ant_time.append(ma.mean(time_limit[0:-2]))
                            ant_mask.append(timebin_mask)
                            ant_limit = []
                            mask_limit = []
                            time_limit = []
                            dataset = 0                            
                        time_limit.append(time)
                        ant_limit.append(new_data)
                        mask_limit.append(new_mask)
                        dataset = dataset + 1
                        if dataset == timescale:
                            timebin_data,timebin_mask = fc.timerebin(ant_limit,mask_limit)
                            rebin_ant.append(timebin_data)
                            ant_time.append(ma.mean(time_limit))
                            ant_mask.append(timebin_mask)
                            ant_limit = []
                            mask_limit = []
                            time_limit = []
                            dataset = 0
                    elif form=='50ohm':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        load.append(new_data)
                        load_time.append(time)
                    elif form=='open':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        term.append(new_data)
                        term_time.append(time)
                    elif form=='short':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        short.append(new_data)
                        short_time.append(time)
                    elif form=='noise':
                        mask = zeros(len(sub_data))
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,30.,200.)
                        noise.append(new_data)
                        noise_time.append(time)
                        
new_freq = array(new_freq)
rebin_ant = array(rebin_ant)
ant_mask= array(ant_mask)
ant_mask_full = zeros((len(ant_mask),len(ant_mask[0])))
for i in range(0,len(new_freq)):
    test_mask = fc.timeflag(rebin_ant[:,i],ant_mask[:,i],ant_time,3.,30)
    ant_mask_full[:,i] = test_mask

median_data = ma.median(rebin_ant,axis=0)

#ant_s11_file = 'Green_Bank_data_may_2013/ANT-120SP-GBT-2013-05-02.s1p'
ant_s11_file = 'Isla_Guadalupe_data_jun_2013/ANT_3_average.s1p'
amp_s_file = 'Isla_Guadalupe_data_jun_2013/WEA101_AMP_2013-04-04.s2p'
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)

data_eff_cal = []
for i in range(0,len(rebin_ant)):
    lin_data = 10.0**(rebin_ant[i]/10.0)
    corr_data = fc.effcal(Eff_sm,new_freq,lin_data)
    data_eff_cal.append(corr_data)

data_eff_cal_db = 10.0*log10(data_eff_cal)
median_eff = ma.median(data_eff_cal_db,axis=0)

R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)
VJF = sqrt(4*1.381e-23*300*50.*(new_freq[1]-new_freq[0])*1e6)*ones(len(new_freq))
VJO = sqrt(4*1.381e-23*300*100.*(new_freq[1]-new_freq[0])*1e6)*ones(len(new_freq))
Z_amp = R_amp_sm(new_freq)+1j*X_amp_sm(new_freq)
Z_ant = R_ant_sm(new_freq)+1j*X_ant_sm(new_freq)
Z50 = 50.*ones(len(new_freq))
Z100 = 100.*exp(2*pi*array(new_freq)*400*1e-6*1j)
#VJAnt = sqrt(4*1.381e-23*300*(new_freq[1]-new_freq[0])*1e6*Z_ant)

In = []
Vn = []
Gain = []
Temp_gain = []
diff_time = []
for j in range(0,len(short_time)):
    for i in range(0,len(load_time)):
        for k in range(0,len(term_time)):
            for l in range(0,len(noise_time)):
                if abs(load_time[i]-short_time[j])<0.01:
                    if abs(load_time[i]-term_time[k])<0.01:
                        if abs(noise_time[l]-load_time[i])<0.01:
                            Vn0,In0,G0,T0 = fc.noise_calc(load[i],short[j],term[k],noise[l],Z_amp,new_freq)
                            Vn.append(Vn0)
                            In.append(In0)
                            Gain.append(G0)
                            Temp_gain.append(T0)
                            diff_time.append(short_time[j])

calscale = 10
calset = 0
In_avg = []
Vn_avg = []
Gain_avg = []
TempG_avg = []
diff_time_avg = []
In_lim = []
Vn_lim = []
Gain_lim = []
TempG_lim = []
diff_time_lim = []

median_Gain = ma.median(Gain,axis=1)
for i in range(0,len(diff_time)):
    if median_Gain[i]>ma.median(median_Gain)/1.25:
        if len(diff_time_lim)<1:
            time_test = diff_time[i]
        else:
            time_test = diff_time_lim[-1]
        if time-time_test>0.5:
            In_avg.append(ma.mean(In_lim,axis=0))
            In_lim = []
            Vn_avg.append(ma.mean(Vn_lim,axis=0))
            Vn_lim = []
            Gain_avg.append(ma.mean(Gain_lim,axis=0))
            Gain_lim = []
            TempG_avg.append(ma.mean(TempG_lim,axis=0))
            TempG_lim = []
            diff_time_avg.append(ma.mean(diff_time_lim,axis=0))
            diff_time_lim = []
            calset=0
        In_lim.append(In[i])
        Vn_lim.append(Vn[i])
        Gain_lim.append(Gain[i])
        TempG_lim.append(Temp_gain[i])
        diff_time_lim.append(diff_time[i])
        calset = calset+1
        if calset==calscale:
            In_avg.append(ma.mean(In_lim,axis=0))
            In_lim = []
            Vn_avg.append(ma.mean(Vn_lim,axis=0))
            Vn_lim = []
            Gain_avg.append(ma.mean(Gain_lim,axis=0))
            Gain_lim = []
            TempG_avg.append(ma.mean(TempG_lim,axis=0))
            TempG_lim = []
            diff_time_avg.append(ma.mean(diff_time_lim,axis=0))
            diff_time_lim = []
            calset=0

if len(diff_time_lim)>=1:
    In_avg.append(ma.mean(In_lim,axis=0))
    In_lim = []
    Vn_avg.append(ma.mean(Vn_lim,axis=0))
    Vn_lim = []
    Gain_avg.append(ma.mean(Gain_lim,axis=0))
    Gain_lim = []
    TempG_avg.append(ma.mean(TempG_lim,axis=0))
    TempG_lim = []
    diff_time_avg.append(ma.mean(diff_time_lim,axis=0))
    diff_time_lim = []
    calset=0    

ant_noise_corr = []
for i in range(0,len(rebin_ant)):
    index1 = 0
    index2 = 0 
    time_comp = 100.
    for j in range(0,len(diff_time)):
        if abs(diff_time[j]-ant_time[i])<time_comp:
            time_comp = abs(diff_time[j]-ant_time[i])
            index2 = index1
            index1 = j
            if abs(diff_time[index2]-ant_time[i])>0.2:
                index2 = j
    if index1==index2
        Psky = fc.noise_corr(rebin_ant[i],Vn[index1],In[index1],new_freq,Z_amp,Z_ant,Gain[index1],Temp_gain[index1])
        Tsky = Temp_gain[index1]*Psky
        ant_noise_corr.append(Tsky)
    else:
        Psky = fc.noise_corr(rebin_ant[i],(Vn[index1]+Vn[index2])/2.,(In[index1]+In[index2])/2.,new_freq,Z_amp,Z_ant,(Gain[index1]+Gain[index2])/2.,(Temp_gain[index1]+Temp_gain[index2])/2.)
        Tsky = (Temp_gain[index1]+Temp_gain[index2])*Psky/2.
        ant_noise_corr.append(Tsky)

median_noise_corr = ma.median(ant_noise_corr,axis=0)

ant_eff_noise_corr = []
#test_eff_noise_corr = []
for i in range(0,len(data_eff_cal_db)):
    index1 = 0
    index2 = 0 
    time_comp = 100.
    for j in range(0,len(diff_time)):
        if abs(diff_time[j]-ant_time[i])<time_comp:
            time_comp = abs(diff_time[j]-ant_time[i])
            index2 = index1
            index1 = j
            if abs(diff_time[index2]-ant_time[i])>0.2:
                index2 = j
    if index1==index2
        Psky = fc.noise_corr(data_eff_cal_db[i],Vn[index1],In[index1],new_freq,Z_amp,Z_ant,Gain[index1],Temp_gain[index1])
        Tsky = Temp_gain[index1]*Psky
        ant_eff_noise_corr.append(Tsky)
    else:
        Psky = fc.noise_corr(data_eff_cal_db[i],(Vn[index1]+Vn[index2])/2.,(In[index1]+In[index2])/2.,new_freq,Z_amp,Z_ant,(Gain[index1]+Gain[index2])/2.,(Temp_gain[index1]+Temp_gain[index2])/2.)
        Tsky = (Temp_gain[index1]+Temp_gain[index2])*Psky/2.
        ant_eff_noise_corr.append(Tsky)

ant_eff_noise_corr = array(ant_eff_noise_corr)
median_eff_noise_corr = ma.median(ant_eff_noise_corr,axis=0)
median_time_eff = ma.median(ant_eff_noise_corr,axis=1)
total_median = ma.median(median_time_eff)

new_mask_full = zeros((len(ant_mask_full),len(ant_mask_full[0])))
for i in range(0,len(ant_mask_full)):
    new_mask_single = fc.flagging(ant_eff_noise_corr[i],new_freq,3.,30)
    new_mask_full[i] = new_mask_single
    spike_mask_single = fc.spike_flag(ant_eff_noise_corr[i],100)
    for j in range(0,len(ant_mask_full[0])):
        if ant_mask_full[i,j]==1.0:
            new_mask_full[i,j] = 1.0
        if median_time_eff[i]>2*total_median:
            new_mask_full[i,j] = 1.0
        if spike_mask_single[j]==1.0:
            new_mask_full[i,j] = 1.0

mean_eff_noise_corr = []
for i in range(0,len(ant_eff_noise_corr[0])):
    test_single = ma.array(ant_eff_noise_corr[:,i],mask=new_mask_full[:,i])
    compress = ma.compressed(test_single)
    mean_single = ma.mean(compress)
    mean_eff_noise_corr.append(mean_single)
    
mean_spike_mask = fc.spike_flag(mean_eff_noise_corr,25)
mean_masked = ma.array(mean_eff_noise_corr,mask = mean_spike_mask)
mean_compressed = ma.compressed(mean_masked)
mean_freq_comp = ma.compressed(ma.array(new_freq,mask=mean_spike_mask))
