"""
Module to process and rebin a single day of data
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
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

##### Set directories
#cal_dir = 'Isla_Guadalupe_data_jun_2013/cal_data/'
cal_dir = '/home/tcv/lustre/cal_data/'
#data_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/June15/'
data_dir = '/home/tcv/lustre/data_arrays/'
out_dir = '/home/tcv/lustre/processed_data_min/'
#out_dir = '/home/tcv/lustre/processed_data/'

##### Load Files
#ant_s11_file = 'Isla_Guadalupe_data_jun_2013/ANT_3_average.s1p'
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
#amp_s_file = 'Isla_Guadalupe_data_jun_2013/WEA101_AMP_2013-04-04.s2p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

diff_time = loadtxt(cal_dir+'Cal_time_take4.txt')
Temp_gain = loadtxt(cal_dir+'TempGain_fit_take4.txt')
GainR = loadtxt(cal_dir+'Real_Gain_fit_take4.txt')
GainX = loadtxt(cal_dir+'Imag_Gain_fit_take4.txt')
InR = loadtxt(cal_dir+'Real_In_fit_take4.txt')
InX = loadtxt(cal_dir+'Imag_In_fit_take4.txt')
VnR = loadtxt(cal_dir+'Real_Vn_fit_take4.txt')
VnX = loadtxt(cal_dir+'Imag_Vn_fit_take4.txt')
new_freq = loadtxt(cal_dir+'Cal_freq_take4.txt')

##### Efficiencey Correction Prep
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)
R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)
Z_amp = R_amp_sm(new_freq)+1j*X_amp_sm(new_freq)
Z_ant = R_ant_sm(new_freq)+1j*X_ant_sm(new_freq)
Z50 = 50.*ones(len(new_freq))
Z100 = 100.*exp(2*pi*array(new_freq)*400*1e-6*1j)


##### Load Data
timescale = 32
binscale = 32
full_rebin_ant = []
full_rebin_time = []
full_rebin_freq = []
full_rebin_volt = []
full_rebin_temp = []

data_arrays = os.listdir(data_dir)
#dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
dates = sys.argv[1]

selected_ind = []
for f in range(0,len(data_arrays)):
    fname = data_arrays[f]
    date = 'June'+dates
    if fname.split('-')[0]==date:
	if fname.split('_')[-1]=='antenna.txt':
            selected_ind.append(f)
	
for s in range(0,len(selected_ind)):
    fname = data_arrays[selected_ind[s]]
#    fname = 'June'+dates[date_ind]+'_antenna.txt'
    if fname.split('_')[-1]=='antenna.txt':
        print 'Data Collection Date is:', fname.split('_')[0]
	sub_date = fname.split('_')[0]
        rebin_ant = loadtxt(data_dir+fname)
        for fname2 in data_arrays:
            if fname2.split('_')[0]==fname.split('_')[0]:
                if fname2.split('_')[-1]=='mask.txt':
                    ant_mask = loadtxt(data_dir+fname2)
                elif fname2.split('_')[-1]=='time.txt':
                    ant_time = loadtxt(data_dir+fname2)
                elif fname2.split('_')[-1]=='freq.txt':
                    new_freq = loadtxt(data_dir+fname2)
                elif fname2.split('_')[-1]=='volt.txt':
                    ant_volt = loadtxt(data_dir+fname2)
                elif fname2.split('_')[-1]=='temp.txt':
                    ant_temp = loadtxt(data_dir+fname2)
        new_freq = array(new_freq)
        rebin_ant = array(rebin_ant)
        ant_mask= array(ant_mask)
        ant_mask_full = zeros((len(ant_mask),len(ant_mask[0])))

##### Flag out Time Dependent Outliers
        for j in range(0,len(new_freq)):
            test_mask = fc.timeflag(rebin_ant[:,j],ant_mask[:,j],ant_time,3.,30)
            ant_mask_full[:,j] = test_mask

##### Correct for Efficiency
        data_eff_cal = []
        for k in range(0,len(rebin_ant)):
            lin_data = 10.0**(rebin_ant[k]/10.0)
            corr_data = fc.effcal(Eff_sm,new_freq,lin_data)
            data_eff_cal.append(corr_data)
        data_eff_cal_db = 10.0*log10(data_eff_cal)
#        data_eff_cal_db = rebin_ant


##### Gain Calibration
        ant_eff_noise_corr = []
        for i in range(0,len(data_eff_cal_db)):
            Tsky = zeros(len(data_eff_cal_db[0]))
            for j in range(0,len(data_eff_cal_db[0])):
               TGain_s = Temp_gain[j]
               Vn_s = VnR[j]+1j*VnX[j]
               In_s = InR[j]+1j*InX[j]
               Gain_s = GainR[j]+1j*GainX[j]
               Tsky[j] = fc.noise_corr(data_eff_cal_db[i,j],Vn_s,In_s,new_freq[j],Z_amp[j],
                                    Z_ant[j],Gain_s,TGain_s)
#               Tsky[j] = TGain_s*Psky
            ant_eff_noise_corr.append(Tsky)
        ant_eff_noise_corr = array(ant_eff_noise_corr)

###### Looking at Time Variation
        single_time_eff = ant_eff_noise_corr[:,4000]
	full_single_median = ma.median(single_time_eff)
        single_time_lim = []
        for i in range(0,len(single_time_eff)):
            if single_time_eff[i]>full_single_median/5.:
		if single_time_eff[i]<5.*full_single_median:
                    single_time_lim.append(single_time_eff[i])            
        median_single_time = ma.median(single_time_lim)
        std_single_time = ma.std(single_time_lim)
	print 'Plotting Time Variation (median is %0.1f)...' %median_single_time 
	
        pylab.scatter(ant_time,single_time_eff,c='b',edgecolor='b',s=3)
        pylab.scatter(ant_time,ones(len(ant_time))*(median_single_time),c='g',edgecolor='g',s=3)
        pylab.scatter(ant_time,ones(len(ant_time))*(median_single_time+3*std_single_time),c='r',edgecolor='r',s=3)
        pylab.scatter(ant_time,ones(len(ant_time))*(median_single_time-3*std_single_time),c='c',edgecolor='c',s=3)
	pylab.legend(('Data','Median','+3 Sigma','-3 Sigma'))
	pylab.grid()
#	pylab.ylim(-2e4,2e4)
	ymin,ymax = pylab.ylim()
	if ymin<0:
	    pylab.ylim(ymin=0)
	if ymax>5e3:
	    pylab.ylim(ymax=5e3)
        pylab.xlim(ant_time[0],ant_time[-1]+0.25)
	pylab.xlabel('Time (Hours since Midnight June 01 GMT)')
	pylab.ylabel('Temperature (Kelvin)')
	pylab.title('Time Variability of Data at %0.1f MHz' %new_freq[4000])
	pylab.savefig(out_dir+sub_date+'_single_freq_median_avgcal',dpi=300)
#        pylab.savefig(data_dir+fname.split('_')[0]+'single_freq_median_avgcal',dpi=300)
        pylab.clf()
   
###### Additional Masking for wild outliers, NaNs etc.
        new_mask_full = zeros((len(ant_mask_full),len(ant_mask_full[0])))
        for i in range(0,len(ant_mask_full)):
            nandata = where(isnan(ant_eff_noise_corr[i]))[0]
            nandata = array(nandata)
            for k in range(0,len(nandata)):
                index = nandata[k]
                new_mask_full[i,index]=1.0
                ant_eff_noise_corr[i,index] = 0.0
            new_mask_single = fc.flagging(ant_eff_noise_corr[i],new_freq,3.,30)
            new_mask_full[i] = new_mask_single
#            spike_mask_single = fc.spike_flag(ant_eff_noise_corr[i],100)
            for j in range(0,len(ant_mask_full[0])):
                if ant_mask_full[i,j]==1.0:
                    new_mask_full[i,j] = 1.0
                if single_time_eff[i]>(median_single_time+3*std_single_time):
                    new_mask_full[i,j] = 1.0
                elif single_time_eff[i]<(median_single_time-3*std_single_time):
                    new_mask_full[i,j] = 1.0
                elif single_time_eff[i]<0:
                    new_mask_full[i,j] = 1.0
#                if single_time_eff[i]<=1000.:
#                    new_mask_full[i,j] = 1.0
#                if spike_mask_single[j]==1.0:
#                    new_mask_full[i,j] = 1.0
	
	print 'Plotting Mask...'
        pylab.imshow(new_mask_full,vmin=-0.1,vmax=1.1,aspect=100./len(new_mask_full),extent=(40,140,len(new_mask_full),0.))
        pylab.ylabel('Time (Sample)')
        pylab.xlabel('Freq (MHz)')
        pylab.colorbar()
        pylab.title('Flagged Frequencies and Times for Data')
	pylab.savefig(out_dir+sub_date+'_data_mask',dpi=300)
#        pylab.savefig(data_dir+fname.split('_')[0]+'_data_mask',dpi=300)
        pylab.clf()
                    
	print 'Plotting Single Time Masking...'
        single_time = int(len(new_mask_full)/2)
	time_used = single_time
	okay = 0
	for j in range(0,single_time):
	    if okay==0:
		print 'Bad File is:',time_used
		meanmask = ma.median(new_mask_full[time_used])
		time_used = time_used-j
		if meanmask<0.9:
		    okay=1
	pylab.scatter(new_freq,ant_eff_noise_corr[time_used],c='b',edgecolor='b',s=3)
	masked_single_time = ma.array(ant_eff_noise_corr[time_used],mask=new_mask_full[time_used])
	compressed_single_time = ma.compressed(masked_single_time)
	masked_single_freq = ma.array(new_freq,mask=new_mask_full[time_used])
	compressed_single_freq = ma.compressed(masked_single_freq)
	print 'Compressed Shape is:',len(compressed_single_time)
	pylab.scatter(compressed_single_freq,compressed_single_time,c='g',edgecolor='g',s=3)
	pylab.grid()
	pylab.xlabel('Frequency (MHz)')
	pylab.xlim(40,140)
	pylab.ylabel('Temperature (Kelvin)')
	pylab.legend(('Unmasked Sample','Masked Sample'))
	pylab.title('Example Masking Dataset')
#	pylab.ylim(-5e4,5e4)
	pylab.ylim(0,2e4)
	pylab.savefig(out_dir+sub_date+'_single_time_masking',dpi=300)
	pylab.clf()


####### Rebin in Frequency 
        freqbin_full = []
        freqbin_mask_full = []
        for i in range(0,len(ant_eff_noise_corr)):
            freqbin_data,freqbin_mask,freqbin_freq = fc.rebin(ant_eff_noise_corr[i],
                                                              new_mask_full[i],
                                                              new_freq,binscale)
            freqbin_full.append(freqbin_data)
            freqbin_mask_full.append(freqbin_mask)

###### Rebin in time
        minimum_time = ant_time[0]
        time_ind = 0
	min_time_ind = 0
        timebin_time = []
	print 'Shape of Time Data is:',shape(ant_time)
        for i in range(0,len(ant_eff_noise_corr)):
	    if time_ind<timescale:
		if time_ind<1:
		    minimum_time = ant_time[i]
		    time_ind = time_ind+1
                elif absolute(ant_time[i]-minimum_time)<1./12:
		    minimum_time = ant_time[i]
                    time_ind = time_ind + 1
		else:
                    print 'Number of Time Samples to Bin:',time_ind
                    timebin_data,timebin_mask = fc.timerebin(freqbin_full[min_time_ind:i],
                                                             freqbin_mask_full[min_time_ind:i])
                    full_rebin_time.append(ma.mean(ant_time[min_time_ind:i]))
                    full_rebin_ant.append(timebin_data)
                    full_rebin_volt.append(ma.mean(ant_volt[min_time_ind:i]))
                    full_rebin_temp.append(ma.mean(ant_temp[min_time_ind:i]))
                    min_time_ind = i+1
                    minimum_time = ant_time[i]
                    time_ind = 0
	    else:
		print 'Number of Time Samples to Bin:',time_ind
                timebin_data,timebin_mask = fc.timerebin(freqbin_full[min_time_ind:i],
                                                         freqbin_mask_full[min_time_ind:i])
                full_rebin_time.append(ma.mean(ant_time[min_time_ind:i]))
                full_rebin_ant.append(timebin_data)
                full_rebin_volt.append(ma.mean(ant_volt[min_time_ind:i]))
                full_rebin_temp.append(ma.mean(ant_temp[min_time_ind:i]))
                min_time_ind = i+1
                minimum_time = ant_time[i]
                time_ind = 0
        if time_ind>=5:
#        if (len(ant_eff_noise_corr)-time_ind)>5:
            timebin_data,timebin_mask = fc.timerebin(freqbin_full[min_time_ind:-1],
                                                     freqbin_mask_full[min_time_ind:-1])
            full_rebin_time.append(ma.mean(ant_time[min_time_ind:-1]))
            full_rebin_ant.append(timebin_data)
            full_rebin_volt.append(ma.mean(ant_volt[min_time_ind:-1]))
            full_rebin_temp.append(ma.mean(ant_temp[min_time_ind:-1]))
        print 'Number of Excess Files:',time_ind
        full_rebin_freq = freqbin_freq
        print 'New Data Array Shape:',shape(full_rebin_ant)

###### Write data to output files
full_rebin_ant = array(full_rebin_ant)
for i in range(0,len(full_rebin_ant)):
    for j in range(0,len(full_rebin_ant[0])):
     	if isnan(full_rebin_ant[i,j])==True:
            full_rebin_ant[i,j]=0.

full_rebin_time = array(full_rebin_time)
full_rebin_freq = array(full_rebin_freq)
full_rebin_volt = array(full_rebin_volt)
full_rebin_temp = array(full_rebin_temp)
savetxt(out_dir+date+'_processed_data_avgcal.txt',full_rebin_ant,delimiter = ' ')
savetxt(out_dir+date+'_processed_time_avgcal.txt',full_rebin_time,delimiter = ' ')
savetxt(out_dir+date+'_processed_freq_avgcal.txt',full_rebin_freq,delimiter = ' ')
savetxt(out_dir+date+'_processed_volt_avgcal.txt',full_rebin_volt,delimiter = ' ')
savetxt(out_dir+date+'_processed_temp_avgcal.txt',full_rebin_temp,delimiter = ' ')
full_rebin_ant = []
full_rebin_time = []
full_rebin_freq = []
full_rebin_volt = []
full_rebin_temp = []
