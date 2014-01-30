import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import sys

#maindir = 'Isla_Guadalupe_data_jun_2013/June15_day_to_night/'
maindir = '/lustre/anat/Guadalupe_data/'
caldir = '/lustre/cal_data_rb/'
#day_dir = 'Isla_Guadalupe_data_jun_2013/cal_data/June15/'
load = []
load_time = []
term = []
term_time = []
short = []
short_time = []
noise = []
noise_time = []
binscale = 1
#timescale = 10

#ant_s11_file = 'Isla_Guadalupe_data_jun_2013/ANT_3_average.s1p'
#amp_s_file = 'Isla_Guadalupe_data_jun_2013/WEA101_AMP_2013-04-04.s2p'
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)

directories = os.listdir(maindir)
date = sys.argv[1]

full_date = 'June'+date
directory = maindir+full_date+'_day_to_night/'
#for direct in directories:
#    full_date = 'June'+date
#    if direct.split('_')[0]==full_date:
#    if direct.split('_')[-1]=='night':
#    if direct.split('-')[0]=='2013':
#        directory = maindir+direct+'/'
#        print shape(load)
for i in range(0,1):
    for j in range(0,1):
	print 'Current Directory Being Processed is:',directory
	load = []
	load_time = []
        load_volt = []
	term = []
	term_time = []
	short = []
	short_time = []
	noise = []
	noise_time = []
	binscale = 1
	day_dir = full_date 
        dirlist = os.listdir(directory)
        for fname in dirlist:
            filename=directory+fname
            if fname.split('.')[-1]=='dat':
                if fname.split('_')[-1]!='antenna.dat':
                    time,form,sub_data,mask,freq,volt = fc.loadsingle(filename)
                    width = 250.0/len(sub_data)
                    freqs = arange(0,250.0,width)
                    if len(sub_data)>1:
                        if form=='50ohm':
                            mask = zeros(len(sub_data))
                            new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                            new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,40.,140.)
                            load.append(new_data)
                            load_time.append(time)
                            load_volt.append(volt)
                        elif form=='open':
                            mask = zeros(len(sub_data))
                            new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                            new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,40.,140.)
                            term.append(new_data)
                            term_time.append(time)
                        elif form=='short':
                            mask = zeros(len(sub_data))
                            new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                            new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,40.,140.)
                            short.append(new_data)
                            short_time.append(time)
                        elif form=='noise':
                            mask = zeros(len(sub_data))
                            new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                            new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,40.,140.)
                            noise.append(new_data)
                            noise_time.append(time)
                        
                        
	new_freq = array(new_freq)
	VJF = sqrt(4*1.381e-23*300*50.*(new_freq[1]-new_freq[0])*1e6)*ones(len(new_freq))
	VJO = sqrt(4*1.381e-23*300*100.*(new_freq[1]-new_freq[0])*1e6)*ones(len(new_freq))
	Z_amp = R_amp_sm(new_freq)+1j*X_amp_sm(new_freq)
	Z_ant = R_ant_sm(new_freq)+1j*X_ant_sm(new_freq)
	Z50 = 50.*ones(len(new_freq))
	Z100 = 100.*exp(2*pi*array(new_freq)*400*1e-6*1j)
        
        gamma_amp = (Z_amp-Z50)/(Z_amp+Z50)
        gamma_ant = (Z_ant-Z50)/(Z_ant+Z50)
        F = sqrt(1-real(gamma_amp*conj(gamma_amp)))/(1-gamma_amp*gamma_ant)
        G = real(1-gamma_amp*conj(gamma_amp))
        gamma_100 = (Z100-Z50)/(Z100+Z50)
        phi = angle(gamma_ant*F)
        F100 = sqrt(1-real(gamma_amp*conj(gamma_amp)))/(1-gamma_amp*gamma_100)
        phi100 = angle(gamma_100*F100)

        Tu = []
        Ts = []
        Tc = []
        int_freq = []
        diff_time = []
        diff_volt = []
        T100_meas = []
	for j in range(0,len(short_time)):
    	    for i in range(0,len(load_time)):
	        for k in range(0,len(term_time)):
        	    for l in range(0,len(noise_time)):
                	if abs(load_time[i]-short_time[j])<0.01:
	                    if abs(load_time[i]-term_time[k])<0.01:
        	                if abs(noise_time[l]-load_time[i])<0.01:
                	            Tu_sin,Ts_sin,Tc_sin,int_freq = fc.rb_noise_calc(load[i],short[j],term[k],noise[l],Z_amp,new_freq)
                        	    diff_time.append(short_time[j])
                                    diff_volt.append(load_volt[i])
                                    Tu.append(Tu_sin)
                                    Ts.append(Ts_sin)
                                    Tc.append(Tc_sin)
                                    T100_sin = (10**(term[i]/10.)-10**(short[j]/10.))/(10**(noise[l]/10.)-10**(load[i]/10.))
                                    T100_meas.append(T100_sin)

	print 'Number of Calibration Datasets is:',shape(diff_time)

        Tu = array(Tu)
	Ts = array(Ts)
	Tc = array(Tc)
        print shape(Tu)

        pylab.plot(new_freq,T100_meas[0])
        pylab.plot(new_freq,T100_meas[1])
        pylab.plot(new_freq,T100_meas[2])
        pylab.plot(new_freq,T100_meas[3])
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Temperature (Kelvin)')
        pylab.title('Test 100 Ohm Resistor Temperature')
        pylab.savefig('/home/tcv/'+caldir+day_dir+'T100_test',dpi=300)
        pylab.clf()
         

        print ma.median(Tu,axis=1)
        print ma.median(Ts,axis=1)
        print ma.median(Tc,axis=1)

	single_Tu = Tu[:,20]
	median_Tu = ma.median(single_Tu)
	std_Tu = ma.std(single_Tu)
        vminTu = 0
        vmaxTu = median_Tu+10*std_Tu
#        if median_Tu<0: 
#            vminTu = median_Tu-10*std_Tu 
#            vmaxTu = median_Tu+10*std_Tu 
#        else: 
#            vminTu = median_Tu-10*std_Tu 
#            vmaxTu = median_Tu+10*std_Tu 
        pylab.imshow(Tu,vmax=vmaxTu,vmin=vminTu,aspect=100./len(Tu),extent=(40,140,len(Tu),0))
        cbar = pylab.colorbar()
        cbar.set_label('Temperature (Noise Cal Temp Units)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')
        pylab.title('Tu Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Tu_waterfall',dpi=300)
        pylab.clf()

        single_Tc = Tc[:,20]
        median_Tc = ma.median(single_Tc)
        std_Tc = ma.std(single_Tc)
        vminTc = 0
        vmaxTc = median_Tc+10*std_Tc
#        if median_Tc<0:
#            vminTc = median_Tc-10*std_Tc
#            vmaxTc = median_Tc+10*std_Tc
#        else:
#            vminTc = median_Tc-10*std_Tc
#            vmaxTc = median_Tc+10*std_Tc
        pylab.imshow(Tc,vmin=vminTc,vmax=vmaxTc,aspect=100./len(Tc),extent=(40,140,len(Tc),0))
        cbar = pylab.colorbar() 
        cbar.set_label('Temperature (Noise Cal Temp Units)') 
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')
        pylab.title('Tc Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Tc_waterfall',dpi=300)
        pylab.clf()

        single_Ts = Ts[:,20]
        median_Ts = ma.median(single_Ts) 
        std_Ts = ma.std(single_Ts) 
        vminTs = 0
        vmaxTs = median_Ts+10*std_Ts
#        if median_Ts<0:
#            vminTs = median_Ts-10*std_Ts
#            vmaxTs = median_Ts+10*std_Ts
#        else:
#            vminTs = median_Ts-10*std_Ts
#            vmaxTs = median_Ts+10*std_Ts

        pylab.imshow(Ts,vmax=vmaxTs,vmin=vminTs,aspect=100./len(Ts),extent=(40,140,len(Ts),0)) 
        cbar = pylab.colorbar()  
        cbar.set_label('Temperature (Noise Cal Temp Units)')  
        pylab.xlabel('Frequency (MHz)') 
        pylab.ylabel('Cal Sample') 
        pylab.title('Ts Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Ts_waterfall',dpi=300)
        pylab.clf() 

	time_sort = argsort(diff_time)
	print time_sort
	Tu_sort = []
	Ts_sort = []
	Tc_sort = []
	diff_time_sort = []
        diff_volt_sort = []
	for t in range(0,len(diff_time)):
	    diff_time_sort.append(diff_time[time_sort[t]])
	    Tu_sort.append(Tu[time_sort[t]])
	    Ts_sort.append(Ts[time_sort[t]])
	    Tc_sort.append(Tc[time_sort[t]])
            diff_volt_sort.append(diff_volt[time_sort[t]])


	savetxt('/home/tcv/'+caldir+day_dir+'Tu_data.txt',Tu_sort,delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Tc_data.txt',Tc_sort,delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Ts_data.txt',Ts_sort,delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Cal_time.txt',diff_time_sort,delimiter=' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Cal_freq.txt',new_freq,delimiter=' ')
        savetxt('/home/tcv/'+caldir+day_dir+'Cal_volt.txt',diff_volt_sort,delimiter=' ')
