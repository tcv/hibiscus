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
caldir = '/lustre/cal_data/'
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
	term = []
	term_time = []
	short = []
	short_time = []
	noise = []
	noise_time = []
	binscale = 1
	day_dir = full_date 
#        if len(direct.split('_'))>2:
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

	print 'Number of Calibration Datasets is:',shape(diff_time)

	if len(Gain)<10.:
	   calscale=len(Gain)
	else:
	   calscale=10
#	calscale = int(len(Gain)/10.)
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

	Gain_array = array(Gain)
	single_Gain = real(Gain_array[:,4000])
	median_Gain = ma.median(real(single_Gain))
	std_Gain = ma.std(real(single_Gain))
	print 'Plotting Gain Evolution over Time...'
	pylab.scatter(diff_time,single_Gain,c='b',edgecolor='b')
	pylab.scatter(diff_time,ones(len(single_Gain))*(median_Gain-3*std_Gain),c='g',edgecolor='g')
	pylab.scatter(diff_time,ones(len(single_Gain))*(median_Gain+3*std_Gain),c='r',edgecolor='r')
#	Gain_fit = itp.UnivariateSpline(diff_time,real(single_Gain))
#	pylab.plot(diff_time,Gain_fit(diff_time))
	pylab.savefig('/home/tcv/'+caldir+day_dir+'_Gain_evolution_take3',dpi=300)
	pylab.clf()

	time_sort = argsort(diff_time)
	print time_sort
	Vn_sort = []
	In_sort = []
	Gain_sort = []
	TempG_sort = []
	diff_time_sort = []
	for t in range(0,len(diff_time)):
	    diff_time_sort.append(diff_time[time_sort[t]])
	    Vn_sort.append(Vn[time_sort[t]])
	    In_sort.append(In[time_sort[t]])
	    Gain_sort.append(Gain[time_sort[t]])
	    TempG_sort.append(Temp_gain[time_sort[t]])

	max_lim = median_Gain+3*std_Gain
	min_lim = median_Gain-3*std_Gain
	if min_lim<0:
	    min_lim = 0

	for i in range(0,len(diff_time_sort)):
	    if single_Gain[i]>min_lim:
        	if len(diff_time_lim)<1:
	            time_test = diff_time_sort[i]
        	else:
	            time_test = diff_time_lim[-1]
#        	if absolute(diff_time[i]-time_test)>0.5:
		if absolute(diff_time_sort[i]-time_test)>1.:
	            if len(In_lim)>1:
                        print 'Number of Samples to be Averaged is:',len(diff_time_lim)
#        	        print 'Number of Averaged Samples is:',len(diff_time_lim)
                	In_lim_mean = ma.mean(In_lim,axis=0)
	                In_avg.append(In_lim_mean)
        	        In_lim = []
                	Vn_lim_mean = ma.mean(Vn_lim,axis=0)
	                Vn_avg.append(Vn_lim_mean)
        	        Vn_lim = []
                	Gain_lim_mean = ma.mean(Gain_lim,axis=0)
	                Gain_avg.append(Gain_lim_mean)
        	        Gain_lim = []
                	TempG_avg.append(ma.mean(TempG_lim,axis=0))
	                TempG_lim = []
        	        diff_time_avg.append(ma.mean(diff_time_lim,axis=0))
                	diff_time_lim = []
	                calset=0
        	In_lim.append(In_sort[i])
	        Vn_lim.append(Vn_sort[i])
        	Gain_lim.append(Gain_sort[i])
	        TempG_lim.append(TempG_sort[i])
        	diff_time_lim.append(diff_time_sort[i])
	        calset = calset+1
        	if calset==calscale:
                    print 'Number of Samples to be Averaged is:',len(diff_time_lim)
#	            print 'Number of Averaged Samples is:',len(diff_time_lim)
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
	    elif single_Gain[i]<max_lim:
        	if len(diff_time_lim)<1:
	            time_test = diff_time_sort[i]
        	else:
	            time_test = diff_time_lim[-1]
#        	if absolute(diff_time[i]-time_test)>0.5:
		if absolute(diff_time_sort[i]-time_test)>1.:
	            if len(In_lim)>1:
                        print 'Number of Samples to be Averaged is:',len(diff_time_lim)
#        	        print 'Number of Averaged Samples is:',len(diff_time_lim)
                	In_lim_mean = ma.mean(In_lim,axis=0)
	                In_avg.append(In_lim_mean)
        	        In_lim = []
                	Vn_lim_mean = ma.mean(Vn_lim,axis=0)
	                Vn_avg.append(Vn_lim_mean)
        	        Vn_lim = []
                	Gain_lim_mean = ma.mean(Gain_lim,axis=0)
	                Gain_avg.append(Gain_lim_mean)
        	        Gain_lim = []
                	TempG_avg.append(ma.mean(TempG_lim,axis=0))
	                TempG_lim = []
        	        diff_time_avg.append(ma.mean(diff_time_lim,axis=0))
                	diff_time_lim = []
	                calset=0
        	In_lim.append(In_sort[i])
	        Vn_lim.append(Vn_sort[i])
        	Gain_lim.append(Gain_sort[i])
	        TempG_lim.append(TempG_sort[i])
        	diff_time_lim.append(diff_time[i])
	        calset = calset+1
        	if calset==calscale:
	            print 'Number of Samples to be Averaged is:',len(diff_time_lim)
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

	if len(diff_time_lim)>1:
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

	print 'Final Number of Averaged Datasets is:',len(diff_time_avg)

	savetxt('/home/tcv/'+caldir+day_dir+'Real_In_avg_take3.txt',real(In_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Imag_In_avg_take3.txt',imag(In_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Real_Vn_avg_take3.txt',real(Vn_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Imag_Vn_avg_take3.txt',imag(Vn_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Real_Gain_avg_take3.txt',real(Gain_avg),delimiter =' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Imag_Gain_avg_take3.txt',imag(Gain_avg),delimiter =' ')
	savetxt('/home/tcv/'+caldir+day_dir+'TempGain_avg_take3.txt',TempG_avg,delimiter=' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Cal_time.txt',diff_time_avg,delimiter=' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Cal_freq.txt',new_freq,delimiter=' ')
