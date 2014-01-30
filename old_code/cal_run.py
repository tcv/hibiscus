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
caldir = '/lustre/cal_data_sept/'
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
        load_temp = []
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
                    time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(filename)
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
                            load_temp.append(temp)
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
        diff_volt = []
        diff_temp = []
	for j in range(0,len(short_time)):
    	    for i in range(0,len(load_time)):
	        for k in range(0,len(term_time)):
        	    for l in range(0,len(noise_time)):
                	if abs(load_time[i]-short_time[j])<0.01:
	                    if abs(load_time[i]-term_time[k])<0.01:
        	                if abs(noise_time[l]-load_time[i])<0.01:
                	            Vn0,In0,G0,T0 = fc.noise_calc(load[i],short[j],term[k],noise[l],Z_amp,new_freq,load_temp[i])
                        	    Vn.append(Vn0)
	                            In.append(In0)
        	                    Gain.append(G0)
                	            Temp_gain.append(T0)
                        	    diff_time.append(short_time[j])
                                    diff_volt.append(load_volt[i])
                                    diff_temp.append(load_temp[i])
   

	print 'Number of Calibration Datasets is:',shape(diff_time)

#if len(Gain)<2.:
#   calscale=len(Gain)
#else:
#   calscale=2
#	calscale = int(len(Gain)/10.)
        calscale=1
	calset = 0
	In_avg = []
	Vn_avg = []
	Gain_avg = []
	TempG_avg = []
	diff_time_avg = []
        diff_volt_avg = []
        diff_temp_avg = []
	In_lim = []
	Vn_lim = []
	Gain_lim = []
	TempG_lim = []
	diff_time_lim = []
        diff_volt_lim = []
        diff_temp_lim = []
        
        Temp_gain = array(Temp_gain)
	print 'Temp Gain Mean is (*1e9):', ma.median(Temp_gain[:,4000])/1e9
        Gain_array = array(Gain)
	single_Gain = real(Gain_array[:,4000])
	median_Gain = ma.median(real(single_Gain))
	std_Gain = ma.std(real(single_Gain))
	print 'Plotting Gain Evolution over Time (median is %0.1f)...' %median_Gain
	pylab.scatter(diff_time,single_Gain,c='b',edgecolor='b')
	pylab.scatter(diff_time,ones(len(single_Gain))*(median_Gain-3*std_Gain),c='g',edgecolor='g')
	pylab.scatter(diff_time,ones(len(single_Gain))*(median_Gain+3*std_Gain),c='r',edgecolor='r')
#	Gain_fit = itp.UnivariateSpline(diff_time,real(single_Gain))
#	pylab.plot(diff_time,Gain_fit(diff_time))
	pylab.savefig('/home/tcv/'+caldir+day_dir+'_Gain_evolution_take4',dpi=300)
	pylab.clf()

        pylab.imshow(real(Gain),vmax=(median_Gain+10*std_Gain),vmin=0,aspect=100./len(Gain),extent=(40,140,len(Gain),0))
        cbar = pylab.colorbar()
        cbar.set_label('Real Gain')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')
        pylab.title('Real Gain Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Real_Gain_waterfall_take4',dpi=300)
        pylab.clf()

        single_GX = imag(Gain_array[:,4000])
        median_GX = ma.median(single_GX)
        std_GX = ma.std(single_GX)

        pylab.imshow(imag(Gain),vmin=(median_GX-10*std_GX),vmax=0,aspect=100./len(Gain),extent=(40,140,len(Gain),0))
        cbar = pylab.colorbar() 
        cbar.set_label('Imaginary Gain') 
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')
        pylab.title('Imag Gain Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Imag_Gain_waterfall_take4',dpi=300)
        pylab.clf()

        single_TG = Temp_gain[:,4000] 
        median_TG = ma.median(single_TG) 
        std_TG = ma.std(single_TG) 

        pylab.imshow(Temp_gain,vmax=(median_TG+100*std_TG),vmin=0,aspect=100./len(Temp_gain),extent=(40,140,len(Temp_gain),0)) 
        cbar = pylab.colorbar()  
        cbar.set_label('Temperature Gain')  
        pylab.xlabel('Frequency (MHz)') 
        pylab.ylabel('Cal Sample') 
        pylab.title('Temperature Gain Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_TempGain_waterfall_take4',dpi=300)
        pylab.clf() 

        Vn = array(Vn)
        single_VR = real(Vn[:,4000])
        median_VR = ma.median(single_VR)
        std_VR = ma.std(single_VR)  
        
        pylab.imshow(real(Vn),vmax=(median_VR+10*std_VR),vmin=0,aspect=100./len(Vn),extent=(40,140,len(Vn),0))
        cbar = pylab.colorbar()   
        cbar.set_label('Real Vn')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')  
        pylab.title('Real Vn Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Real_Vn_waterfall_take4',dpi=300)
        pylab.clf()  

        single_VX = imag(Vn[:,4000])
        median_VX = ma.median(single_VX)
        std_VX = ma.std(single_VX)  
        
        pylab.imshow(imag(Vn),vmax=(median_VX+10*std_VX),vmin=0,aspect=100./len(Vn),extent=(40,140,len(Vn),0))
        cbar = pylab.colorbar()   
        cbar.set_label('Imag Vn')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')  
        pylab.title('Imag Vn Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Imag_Vn_waterfall_take4',dpi=300)
        pylab.clf()  

        In = array(In)
        single_IR = real(In[:,4000])
        median_IR = ma.median(single_IR)
        std_IR = ma.std(single_IR)  
        
        pylab.imshow(real(In),vmax=(median_IR+10*std_IR),vmin=0,aspect=100./len(In),extent=(40,140,len(In),0))
        cbar = pylab.colorbar()   
        cbar.set_label('Real In')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')  
        pylab.title('Real In Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Real_In_waterfall_take4',dpi=300)
        pylab.clf()  

        single_IX = imag(In[:,4000])
        median_IX = ma.median(single_IX)
        std_IX = ma.std(single_IX)  
        
        pylab.imshow(imag(In),vmax=(median_IX+10*std_IX),vmin=0,aspect=100./len(In),extent=(40,140,len(In),0))
        cbar = pylab.colorbar()   
        cbar.set_label('Imag In')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Cal Sample')  
        pylab.title('Imaginary In Data Variation over time for '+ full_date)
        pylab.savefig('/home/tcv/'+caldir+day_dir+'_Imag_In_waterfall_take4',dpi=300)
        pylab.clf()  

	time_sort = argsort(diff_time)
	print time_sort
	Vn_sort = []
	In_sort = []
	Gain_sort = []
	TempG_sort = []
	diff_time_sort = []
        diff_volt_sort = []
        diff_temp_sort = []
	for t in range(0,len(diff_time)):
	    diff_time_sort.append(diff_time[time_sort[t]])
	    Vn_sort.append(Vn[time_sort[t]])
	    In_sort.append(In[time_sort[t]])
	    Gain_sort.append(Gain[time_sort[t]])
	    TempG_sort.append(Temp_gain[time_sort[t]])
            diff_volt_sort.append(diff_volt[time_sort[t]])
            diff_temp_sort.append(diff_temp[time_sort[t]])

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
	            if len(In_lim)>=1:
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
                        diff_volt_avg.append(ma.mean(diff_volt_lim,axis=0))
                        diff_volt_lim = []
                        diff_temp_avg.append(ma.mean(diff_temp_lim,axis=0))
                        diff_temp_lim = []
	                calset=0
        	In_lim.append(In_sort[i])
	        Vn_lim.append(Vn_sort[i])
        	Gain_lim.append(Gain_sort[i])
	        TempG_lim.append(TempG_sort[i])
        	diff_time_lim.append(diff_time_sort[i])
                diff_volt_lim.append(diff_volt_sort[i])
                diff_temp_lim.append(diff_temp_sort[i])
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
                    diff_volt_avg.append(ma.mean(diff_volt_lim,axis=0))
                    diff_volt_lim = []
                    diff_temp_avg.append(ma.mean(diff_temp_lim,axis=0))
                    diff_temp_lim = []
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
                        diff_volt_avg.append(ma.mean(diff_volt_lim,axis=0))
                        diff_volt_lim = []
                        diff_temp_avg.append(ma.mean(diff_temp_lim,axis=0))
                        diff_temp_lim = []
	                calset=0
        	In_lim.append(In_sort[i])
	        Vn_lim.append(Vn_sort[i])
        	Gain_lim.append(Gain_sort[i])
	        TempG_lim.append(TempG_sort[i])
        	diff_time_lim.append(diff_time_sort[i])
                diff_volt_lim.append(diff_volt_sort[i])
                diff_temp_lim.append(diff_temp_sort[i])
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
                    diff_volt_avg.append(ma.mean(diff_volt_lim,axis=0))
                    diff_volt_lim = []
                    diff_temp_avg.append(ma.mean(diff_temp_lim,axis=0))
                    diff_temp_lim = []
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
            diff_volt_avg.append(ma.mean(diff_volt_lim,axis=0))
            diff_volt_lim = []
            diff_temp_avg.append(ma.mean(diff_temp_lim,axis=0))
            diff_temp_lim = []
	    calset=0

	print 'Final Number of Averaged Datasets is:',len(diff_time_avg)

	savetxt('/home/tcv/'+caldir+day_dir+'Real_In.txt',real(In_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Imag_In.txt',imag(In_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Real_Vn.txt',real(Vn_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Imag_Vn.txt',imag(Vn_avg),delimiter = ' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Real_Gain.txt',real(Gain_avg),delimiter =' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Imag_Gain.txt',imag(Gain_avg),delimiter =' ')
	savetxt('/home/tcv/'+caldir+day_dir+'TempGain.txt',TempG_avg,delimiter=' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Cal_time.txt',diff_time_avg,delimiter=' ')
	savetxt('/home/tcv/'+caldir+day_dir+'Cal_freq.txt',new_freq,delimiter=' ')
        savetxt('/home/tcv/'+caldir+day_dir+'Cal_volt.txt',diff_volt_avg,delimiter=' ')
        savetxt('/home/tcv/'+caldir+day_dir+'Cal_temp.txt',diff_temp_avg,delimiter=' ')
