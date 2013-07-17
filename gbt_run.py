from numpy import *
import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf

olddir = 'Green_Bank_data_may_2013/Redwood_house_data/'
old_rebin_ant =[]
old_ant_time = []
old_mask = []
old_load = []
old_load_time = []
old_term = []
old_term_time = []
old_short = []
old_short_time = []
old_noise = []
old_noise_time = []
old_directories = os.listdir(olddir)
for direct in old_directories:
    if direct.split('-')[0]=='2013':
        directory = olddir+direct+'/'
        print directory
        print shape(old_rebin_ant)
        dirlist = os.listdir(directory)
        for fname in dirlist:
            filename=directory+fname
            time,form,sub_data,mask,freq,volt = fc.loadsingle(filename)
            width = 250.0/len(sub_data)
            freqs = arange(0,250.0,width)
            binscale=25
            if len(sub_data)>1:
                if form=='antenna':
#                mask = zeros(len(sub_data))
                    mask = fc.flagging(sub_data,freqs,3.,1e4)
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    old_rebin_ant.append(new_data)
                    old_ant_time.append(time)
                    old_mask.append(new_mask)
                elif form=='50ohm':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    old_load.append(new_data)
                    old_load_time.append(time)
                elif form=='open':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    old_term.append(new_data)
                    old_term_time.append(time)
                elif form=='short':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    old_short.append(new_data)
                    old_short_time.append(time)
                elif form=='noise':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    old_noise.append(new_data)
                    old_noise_time.append(time)
        
old_rebin_freqs = new_freq

maindir = 'Green_Bank_data_may_2013/Old_antenna_house_data/'
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

directories = os.listdir(maindir)
for direct in directories:
    if direct.split('-')[0]=='2013':
        directory = maindir+direct+'/'
        print directory
        print shape(rebin_ant)
        dirlist = os.listdir(directory)
        for fname in dirlist:
            filename=directory+fname
            time,form,sub_data,mask,freq,volt = fc.loadsingle(filename)
            width = 250.0/len(sub_data)
            freqs = arange(0,250.0,width)
            binscale=25
            if len(sub_data)>1:
                if form=='antenna':
#                mask = zeros(len(sub_data))
                    mask = fc.flagging(sub_data,freqs,3.,1e4)
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    rebin_ant.append(new_data)
                    ant_time.append(time)
                    ant_mask.append(new_mask)
                elif form=='50ohm':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    load.append(new_data)
                    load_time.append(time)
                elif form=='open':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    term.append(new_data)
                    term_time.append(time)
                elif form=='short':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    short.append(new_data)
                    short_time.append(time)
                elif form=='noise':
                    mask = zeros(len(sub_data))
                    new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                    noise.append(new_data)
                    noise_time.append(time)
                

median_data = ma.median(rebin_ant,axis=0)
old_median_data = ma.median(old_rebin_ant,axis=0)

ant_s11_file = 'Green_Bank_data_may_2013/ANT-120SP-GBT-2013-05-02.s1p'
amp_s_file = 'Green_Bank_data_may_2013/WEA101_AMP_2013-04-04.s2p'
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)

Eff_sm = fc.smooth(Effic,F_ant,0.01)

old_data_effcal = []
data_eff_cal = []
for i in range(0,len(rebin_ant)):
    lin_data = 10.0**(rebin_ant[i]/10.0)
    corr_data = fc.effcal(Eff_sm,new_freq,lin_data)
    data_eff_cal.append(corr_data)

ant_s11_file_old = 'Green_Bank_data_may_2013/ANT-120SP-GBT-2013-05-01.s1p'
R_ant_old,X_ant_old,F_ant_old = fc.imped_skrf(ant_s11_file_old,0.0)
Eff_old = fc.effic(R_ant_old,X_ant_old,F_ant_old,R_amp,X_amp,F_amp)
Eff_sm_old = fc.smooth(Eff_old,F_ant_old,0.01)

for i in range(0,len(old_rebin_ant)):
    old_lin_data = 10.0**(old_rebin_ant[i]/10.0)
    old_corr_data = fc.effcal(Eff_sm_old,old_rebin_freqs,old_lin_data)
    old_data_effcal.append(old_corr_data)

data_eff_cal_db = 10.0*log10(data_eff_cal)
old_data_effcal_db = 10*log10(old_data_effcal)
median_eff = ma.median(data_eff_cal_db,axis=0)
median_eff_old=ma.median(old_data_effcal_db,axis=0)

cal_list = os.listdir('Green_Bank_data_may_2013/Cal_data/')
cal_data = []
cal_time = []
cal_form = []
for fname in cal_list:
    filename = 'Green_Bank_data_may_2013/Cal_data/'+fname
    time,form,data,mask,freq,volt = fc.loadsingle(filename)
    mask = zeros(len(data))
    new_data,new_mask,new_freq = fc.rebin(data,mask,freqs,binscale)
    cal_data.append(new_data)
    cal_time.append(time)
    cal_form.append(form)

R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)
VJF = sqrt(4*1.381e-23*300*50.*(new_freq[1]-new_freq[0])*1e6)*ones(len(new_freq))
VJO = sqrt(4*1.381e-23*300*100.*(new_freq[1]-new_freq[0])*1e6)*ones(len(new_freq))
Z_amp = R_amp_sm(new_freq)+1j*X_amp_sm(new_freq)
Z_ant = R_ant_sm(new_freq)+1j*X_ant_sm(new_freq)
Z50 = 50.*ones(len(new_freq))
Z100 = 100.*e**(2*pi*new_freq*400*1e-6*1j)
#VJAnt = sqrt(4*1.381e-23*300*(new_freq[1]-new_freq[0])*1e6*Z_ant)
#NF = 5-0*log10(10**((cal_data[1]-cal_data[2])/10.)-1)
#NF_red = 5-10*log10(10**((old_noise[7]-old_load[7])/10.)-1)

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

Pnoise = 10**(cal_data[1]/10.)
P50 = 10**(cal_data[2]/10.)
Tns = 300*(Pnoise/P50)
Vns = (Z_amp+Z50)*sqrt(Z_amp)*(sqrt(Pnoise)-sqrt(P50))/Z_amp
Pns = absolute(Vns**2/Z_amp)
#Pns = (Pnoise-P50)
Pnoise_red = 10**(old_noise[7]/10.)
P50_red = 10**(old_load[7]/10.)
Tns_red = 300*(Pnoise_red/P50_red)
Vns_red = (Z_amp+Z50)*sqrt(Z_amp)*(sqrt(Pnoise_red)-sqrt(P50_red))/Z_amp
Pns_red = absolute(Vns_red**2/Z_amp)
#Pns_red = Pnoise_red-P50_red

ant_noise_corr = []
ant_noise_corr_red = []
for i in range(0,len(rebin_ant)):
    index = 0
    time_comp = 100.
    for j in range(0,len(diff_time)):
        if abs(diff_time[j]-ant_time[i])<time_comp:
            time_comp = abs(diff_time[j]-ant_time[i])
            index = j
    Psky = fc.noise_corr(rebin_ant[i],Vn[index],In[index],new_freq,Z_amp,Z_ant,Gain[index],Temp_gain[index])
    Tsky = Psky*Tns/(Pns/absolute(Gain[index]**2))
    Tsky_red = Psky*Tns_red/(Pns_red/absolute(Gain[index]**2))
    ant_noise_corr.append(Tsky)
    ant_noise_corr_red.append(Tsky_red)

median_noise_corr = ma.median(ant_noise_corr,axis=0)
median_noise_corr_red = ma.median(ant_noise_corr_red,axis=0)

ant_eff_noise_corr = []
ant_eff_noise_corr_red = []
for i in range(0,len(data_eff_cal_db)):
    index = 0
    time_comp = 100.
    for j in range(0,len(diff_time)):
        if abs(diff_time[j]-ant_time[i])<time_comp:
            time_comp = abs(diff_time[j]-ant_time[i])
            index = j
    Psky = fc.noise_corr(data_eff_cal_db[i],Vn[index],In[index],new_freq,Z_amp,Z_ant,Gain[index],Temp_gain[index])
    Tsky = Psky*Tns/(Pns/absolute(Gain[index]**2))
    Tsky_red = Psky*Tns_red/(Pns_red/absolute(Gain[index]**2))
    ant_eff_noise_corr.append(Tsky)
    ant_eff_noise_corr_red.append(Tsky_red)

median_eff_noise_corr = ma.median(ant_eff_noise_corr,axis=0)
median_eff_noise_corr_red = ma.median(ant_eff_noise_corr_red,axis=0)

#ant_noise_corr_redwood = []
#for i in range(0,len(rebin_ant)):
#    index = 0
#    time_comp = 100.
#    for j in range(0,len(diff_time)):
#        if abs(diff_time[j]-ant_time[i])<time_comp:
#            time_comp = abs(diff_time[j]-ant_time[i])
#            index = j
#    Psky = fc.noise_corr(rebin_ant[i],Vn_red[index],In50_red[index],new_freq,Z_amp,Z_ant,Gain_diff_redwood_sm(new_freq))
#    Tsky = Psky*gain_kelvin_red
#    ant_noise_corr_redwood.append(Tsky)

#median_noise_corr_redwood = ma.median(ant_noise_corr_redwood,axis=0)

#ant_eff_noise_redwood = []
#for i in range(0,len(data_eff_cal_db)):
#    index = 0
#    time_comp = 100.
#    for j in range(0,len(diff_time)):
#        if abs(diff_time[j]-ant_time[i])<time_comp:
#            time_comp = abs(diff_time[j]-ant_time[i])
#            index = j
#    Psky = fc.noise_corr(data_eff_cal_db[i],Vn_red[index],In50_red[index],new_freq,Z_amp,Z_ant,Gain_diff_redwood_sm(new_freq))
#    Tsky = Psky*gain_kelvin_red
#    ant_eff_noise_redwood.append(Tsky)

#median_eff_noise_redwood = ma.median(ant_eff_noise_redwood,axis=0)
#Gain = []
#Gain_time = []
#for i in range(0,len(term_time)):
#    index=0
#    time_comp=100.
#    for j in range(0,len(diff_time)):
#        if abs(diff_time[j]-term_time[i])<time_comp:
#            time_comp=abs(diff_time[j]-term_time[i])
#            index=j
#    G0 = fc.gain_calc_term(term[i],Vn[index],In50[index],Z_amp,new_freq,Z100)
#    Gain.append(G0)
#    Gain_time.append(term_time[i])

#Gain_corr_term =[]
#for i in range(0,len(ant_noise_rm)):
#    index = 0
#    time_comp = 100.
#    for j in range(0,len(Gain_time)):
#        if abs(Gain_time[j]-ant_time[i])<time_comp:
#            time_comp = abs(Gain_time[j]-ant_time[i])
#            index = j
#    corr_data = fc.gain_corr_term(Gain[index],ant_noise_rm[i])
#    Gain_corr_term.append(corr_data)

#median_gainterm = ma.median(Gain_corr_term,axis=0)

#R_ant_old_sm = fc.smooth(R_ant_old,F_ant_old,0)
#X_ant_old_sm = fc.smooth(X_ant_old,F_ant_old,0)
#Z_ant_old = R_ant_old_sm(new_freq)+1j*X_ant_old_sm(new_freq)
#VJAnt_old = sqrt(4*1.381e-23*300*(new_freq[1]-new_freq[0])*1e6*Z_ant_old)

#In50_old = []
#In100_old = []
#Vn_old = []
#diff_time_old = []
#for j in range(0,len(old_short_time)):
#    for i in range(0,len(old_load_time)):
#        if abs(old_load_time[i]-old_short_time[j])<0.005:
#            Vn0,In0=fc.noise_calc(old_load[i],old_short[j],Z_amp,new_freq,Z50)
#            Vn_old.append(Vn0)
#            In50_old.append(In0)
#            diff_time_old.append(old_short_time[j])
#    for i in range(0,len(old_term_time)):
#        if abs(old_term_time[i]-old_short_time[j])<0.005:
#            Vn0,In0=fc.noise_calc(old_term[i],old_short[j],Z_amp,new_freq,Z100)
#            In100_old.append(In0)

#old_ant_noise_rm = []
#for i in range(0,len(old_rebin_ant)):
#    index = 0
#    time_comp = 100.
#    for j in range(0,len(diff_time_old)):
#        if abs(diff_time_old[j]-old_ant_time[i])<time_comp:
#            time_comp = abs(diff_time_old[j]-old_ant_time[i])
#            index = j
#    Psky_db = fc.noise_corr(old_rebin_ant[i],Vn_old[index],In50_old[index],new_freq,Z_amp,Z_ant_old)
#    old_ant_noise_rm.append(Psky_db)

#old_Gain = []
#old_Gain_time = []
#for i in range(0,len(old_term_time)):
#    index=0
#    time_comp=100.
#    for j in range(0,len(diff_time_old)):
#        if abs(diff_time_old[j]-old_term_time[i])<time_comp:
#            time_comp=abs(diff_time_old[j]-old_term_time[i])
#            index=j
#    G0 = fc.gain_calc_term(old_term[i],Vn_old[index],In50_old[index],Z_amp,new_freq,Z100)
#    old_Gain.append(G0)
#    old_Gain_time.append(old_term_time[i])

#old_Gain_corr_term =[]
#for i in range(0,len(old_ant_noise_rm)):
#    index = 0
#    time_comp = 100.
#    for j in range(0,len(old_Gain_time)):
#        if abs(old_Gain_time[j]-old_ant_time[i])<time_comp:
#            time_comp = abs(old_Gain_time[j]-old_ant_time[i])
#            index = j
#    corr_data = fc.gain_corr_term(old_Gain[index],old_ant_noise_rm[i])
#    old_Gain_corr_term.append(corr_data)

#old_median_gainterm = ma.median(old_Gain_corr_term,axis=0)
