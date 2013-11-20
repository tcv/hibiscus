"""
Module for plotting the processed data outputs.
"""

import matplotlib
matplotlib.use('Agg')
from numpy import *
import matplotlib.pyplot as pylab
#from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
from scipy import linalg
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import ephem as eph

#Initial Settings:
gsm = False #If want to do GSM model comparison
full = False #If want to expand waterfall to fill in gaps
eff = True #If want to play with efficiency and time delay to uneff data
fitting = True #If want to do power law fitting of data means


#Read in processed data
#data_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/'
data_dir = '/home/tcv/lustre/processed_data_sept/'
#result_dir = '/home/tcv/lustre/data_plots_sept/'
result_dir = '/home/tcv/data_plots_sept_take2/'
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_guad_take2.dat')
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

data = os.listdir(data_dir)

processed_data = []
processed_time = []
processed_freq = []
processed_volt = []
dates = ['01','03','04','05']
#dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
date_files = []
day_ends = []
day_ends.append(0)

if gsm:
    gsm_freq = arange(60,90,1)
    gsm_time = arange(0,24,24./289)
    gsm_data = zeros((len(gsm_time),len(gsm_freq)))
    for t in range(0,len(gsm_time)):
        for f in range(0,len(gsm_freq)):
            gsm_data[t,f] = gsm_raw_data[f*289+t,2]

for i in range(0,len(dates)):
    curr_date = 'June'+dates[i]
    print 'Date being processed is:',curr_date
    single_data = loadtxt(data_dir+curr_date+'_processed_data_subavg.txt')
    single_data = array(single_data)
    for k in range(0,len(single_data)):
	processed_data.append(single_data[k])
    single_time = loadtxt(data_dir+curr_date+'_processed_time_subavg.txt')
    for j in range(0,len(single_time)):
        processed_time.append(single_time[j])
    single_volt = loadtxt(data_dir+curr_date+'_processed_volt_subavg.txt')
    for l in range(0,len(single_volt)):
        processed_volt.append(single_volt[l])
    if len(processed_freq)<1:
        processed_freq = loadtxt(data_dir+curr_date+'_processed_freq_subavg.txt')
    print 'Current Size of Processed Data is:',shape(processed_data)    
    day_ends.append(len(processed_data))
	    
processed_data = array(processed_data)
processed_time = array(processed_time)
processed_volt = array(processed_volt)
print 'Full Shape of Processed Data is:',shape(processed_data)
print 'Shape of Processed Time is:',shape(processed_time)
print 'Nan Data present',where(isnan(processed_data))[0]
print 'Inf Data Present',where(isinf(processed_data))[0]

#Limited Frequencies to Try
#data_lim = [55.,90.]
data_lim = [60.,90.]
##################################################################################
data_lim_ind = [where(processed_freq<=data_lim[0])[0][-1],where(processed_freq>=data_lim[1])[0][0]]

# Adding Masking
prelim_mask = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(prelim_mask)):
    spike_mask = fc.spike_flag(processed_data[i],100)
    for j in range(0,len(prelim_mask[0])):
	if spike_mask[j] == 1.0:
	    prelim_mask[i,j] = 1.0
            processed_data[i,j] = ma.mean(processed_data[i,j-5:j+5])
        if processed_data[i,j]==0.0:
            prelim_mask[i,j] = 1.0
            processed_data[i,j] = ma.mean(processed_data[i,j-5:j+5])
for i in range(0,len(prelim_mask[0])):
    spike_mask = fc.spike_flag(processed_data[:,i],100)
    for j in range(0,len(prelim_mask)):
        if spike_mask[j]==1.0:
            prelim_mask[j,i] = 1.0
            processed_data[j,i] = ma.mean(processed_data[j-5:j+5,i])
for i in range(0,len(prelim_mask)):
    if processed_volt[i]<=9.5:
        for j in range(0,len(processed_freq)):
            prelim_mask[i,j]==1.0
    elif processed_volt[i]>=12.0:
        for j in range(0,len(processed_freq)):
            prelim_mask[i,j]==1.0

#pylab.imshow(prelim_mask,aspect=1.0*len(prelim_mask[0])/len(prelim_mask))
#pylab.xlabel('Frequency Bin')
#pylab.ylabel('Time Bin')
#pylab.title('Data Mask')
#pylab.savefig(result_dir+'prelim_mask',dpi=300)
#pylab.clf()

battery_cuts = []
battery_cuts.append(0)
for i in range(1,len(processed_volt)):
    if (processed_volt[i]-processed_volt[i-1])>1.0:
        battery_cuts.append(i)

print 'Cuts for Battery Swapping:',battery_cuts

#IF need to add efficiency data to other data
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic)

#Base level Mean Data
mean_data = []
for i in range(0,len(processed_data[0])):
    single_freq = ma.array(processed_data[:,i],mask=prelim_mask[:,i])
    single_compress = ma.compressed(single_freq)
    single_mean = ma.mean(single_compress)
    mean_data.append(single_mean)
nandata = where(isnan(mean_data))
nandata = array(nandata)    
for i in range(0,len(nandata[0])):
    index=nandata[0,i]
    mean_data[index]=0.01
infdata = where(isinf(mean_data))
infdata = array(infdata)
for i in range(0,len(infdata[0])):
    index=infdata[0,i]
    mean_data[index]=0.01

#Battery level Mean Data
battery_mean_data = zeros((len(battery_cuts)-1,len(processed_freq)))
for i in range(1,len(battery_cuts)):
    for j in range(0,len(processed_freq)):
        single_freq = ma.array(processed_data[battery_cuts[i-1]:battery_cuts[i],j],
                               mask=prelim_mask[battery_cuts[i-1]:battery_cuts[i],j])
        single_compress = ma.compressed(single_freq)
        single_mean= ma.mean(single_compress)
        battery_mean_data[i-1,j]=single_mean
    nandata = where(isnan(battery_mean_data[i-1]))
    nandata = array(nandata)    
    for i in range(0,len(nandata[0])):
        index=nandata[0,i]
        battery_mean_data[i-1,index]=0.01
    infdata = where(isinf(battery_mean_data[i-1]))
    infdata = array(infdata)
    for i in range(0,len(infdata[0])):
        index=infdata[0,i]
        battery_mean_data[i-1,index]=0.01
    

#zerodata = where(mean_data==0.0)
#zerodata=array(zerodata)
#for i in range(0,len(zerodata[0])):
#    index=zerodata[0,i]
#    mean_data[index] = 0.01


#lgdata = where(mean_data>2e4)
#lgdata = array(lgdata)
#for i in range(0,len(lgdata[0])):
#    index = lgdata[0,i]

#Adding efficiency component
if eff:
    time_delay = arange(0,1.4e-9,2.e-10)
    Eff_sm_mean = []
    R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
    for i in range(0,len(time_delay)):
        R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,time_delay[i])
        Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
        Eff_sm = itp.UnivariateSpline(F_ant,Effic)
        new_mean_data = fc.effcal(Eff_sm,processed_freq,mean_data)
        Eff_sm_mean.append(new_mean_data)
    Eff_sm_mean = array(Eff_sm_mean)
    colors = ['b','g','r','c','m','y','k']
    for i in range(0,len(time_delay)):
	pylab.plot(processed_freq,Eff_sm_mean[i],c=colors[i],label=str(time_delay[i]*1.e9))
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Temperature (Kelvin)')
    pylab.title('Efficiency Corrected Data')
    pylab.ylim(0,2e4)
    pylab.legend()
    pylab.grid()
    pylab.savefig(result_dir+'mean_freq_spectrum_with_eff_compare_june1-5',dpi=300)
    pylab.clf()

    if fitting:
        fitfunceff = lambda x,p0,p1,p2: p0+p1*log10(x)+p2*(log10(x))**2
        piniteff = [-2.5,11,0.1]
        eff_mean_fits = zeros((len(time_delay),3))
        eff_resid = []
        for i in range(0,len(time_delay)):
            peff,success = opt.curve_fit(fitfunceff,processed_freq[data_lim_ind[0]:data_lim_ind[1]],log10(Eff_sm_mean[i,data_lim_ind[0]:data_lim_ind[1]]),piniteff[:],maxfev=1000000)
            eff_mean_fits[i] = peff
            eff_resid.append(Eff_sm_mean[i] - 10**(fitfunceff(processed_freq,peff[0],peff[1],peff[2])))
        eff_resid = array(eff_resid)
        eff_mean_fits = array(eff_mean_fits)
        print 'Residual Shape is:',shape(eff_resid)
        
        lim_eff_resid = eff_resid[:,data_lim_ind[0]:data_lim_ind[1]]
        mean_resid = []
        for i in range(0,len(lim_eff_resid)):
            mean_single = ma.mean(absolute(lim_eff_resid[i]))
            mean_resid.append(mean_single)
        print 'Number of Residual datasets is:',len(mean_resid)
        print 'Time Delay     Mean Residual'
        for i in range(0,len(mean_resid)):
            print time_delay[i],mean_resid[i]
        min_time_delay = where(mean_resid==amin(mean_resid))[0]
        print 'Minimum Residual is at:', time_delay[min_time_delay[0]]
        print 'Minimum Residual Parameters are:', eff_mean_fits[min_time_delay[0]]

        pylab.scatter(processed_freq,Eff_sm_mean[min_time_delay[0]],c='b',edgecolor='b',s=1,label='Eff Corrected Data')
        pylab.plot(processed_freq,10**(fitfunceff(processed_freq,eff_mean_fits[min_time_delay[0],0],eff_mean_fits[min_time_delay[0],1], eff_mean_fits[min_time_delay[0],2])),c='g',label='Fit')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Temperature (Kelvin)')
        pylab.ylim(0,2e4)
        pylab.grid()
        delay = float(time_delay[min_time_delay[0]])*1e9
        pylab.title('Best Fit using Efficiency with Phase Delay')
#        of %0.1f ns' %delay)
        pylab.savefig(result_dir +'mean_freq_spectrum_with_opteff_fit_june1-5',dpi=300)
        pylab.clf()

        pylab.scatter(processed_freq,eff_resid[min_time_delay[0]],c='b',edgecolor='b',s=1,label='best resid, %0.1f ns' %time_delay[min_time_delay[0]])
        pylab.ylim(-500,500)
        pylab.grid()
        pylab.ylabel('Temperature (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.title('Residuals for 3 param power law fit to Efficiency Corrected Data at %0.1f-%0.1f MHz' %(data_lim[0],data_lim[1]))
        pylab.savefig(result_dir+'mean_freq_spectrum_with_opteff_resid_june1-5',dpi=300)
        pylab.clf()
          
#Base Level power law fit
if fitting:
    fitfunc = lambda x,p0,p1,p2: p0+p1*log10(x)+p2*log10(x)**2
#rrfunc = lambda p,x,y: fitfunc(p,x)-y
    pinit = [11,-2.5,1.]
    p,success = opt.curve_fit(fitfunc,processed_freq,log10(mean_data),pinit[:],ones(len(mean_data)),maxfev=1000000)
    pinit = p
    p,success = opt.curve_fit(fitfunc,processed_freq,log10(mean_data),pinit[:],ones(len(mean_data)),maxfev=1000000)
#p,success = opt.leastsq(errfunc,p0,args=(processed_freq,mean_data))
    p1,success1 = opt.curve_fit(fitfunc,processed_freq[data_lim_ind[0]:data_lim_ind[1]],log10(mean_data[data_lim_ind[0]:data_lim_ind[1]]),pinit[:],ones(len(mean_data[data_lim_ind[0]:data_lim_ind[1]])),maxfev=10000000)
    pinit = p1
#p1,success1 = opt.leastsq(errfunc,p0,args=(processed_freq[data_lim_ind[0]:data_lim_ind[1]],mean_data[data_lim_ind[0]:data_lim_ind[1]]))
    fitfunc2 = lambda x,q0,q1: q0-2.5*log10(x)**(-2.5)+q1*log10(x)**2
#errfunc2 = lambda q,a,b: fitfunc2(q,a)-b
    qinit = [11.,1.]
    q,success2 = opt.curve_fit(fitfunc2,processed_freq,log10(mean_data),qinit[:],ones(len(mean_data)))
    q1,success3= opt.curve_fit(fitfunc2,processed_freq[data_lim_ind[0]:data_lim_ind[1]],log10(mean_data[data_lim_ind[0]:data_lim_ind[1]]),qinit[:],ones(len(mean_data[data_lim_ind[0]:data_lim_ind[1]])))
#q,success2 = opt.leastsq(errfunc2,q0,args=(processed_freq,mean_data))
#q1,success3 = opt.leastsq(errfunc2,q0,args=(processed_freq[data_lim_ind[0]:data_lim_ind[1]],mean_data[data_lim_ind[0]:data_lim_ind[1]]))

    print 'Full Data Fit Params (variable power law) are:',p
    print 'Full Data Fit Param is:',q
    print 'Limited Freq Data Fit Params (variable power law) are:',p1
    print 'Limited Freq Data Fit Param is:', q1


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

if full:
    full_time = arange(processed_time[0],processed_time[-1],0.01)
    full_data = zeros((len(full_time),len(processed_freq)))
    for i in range(0,len(full_data)):
        for j in range(0,len(processed_time)):
            if abs(full_time[i]-processed_time[j])<=0.01:
                full_data[i]=processed_data[j]
    
    full_time_sidereal = []
    for t in range(0,len(full_time)):
        single_date = eph.date(initial+full_time[t]/24.)
        guad.date = single_date
        single_time = guad.sidereal_time()
        single_hour = single_time*12./pi
        full_time_sidereal.append(single_hour)

    if gsm:
        full_gsm_data = zeros((len(full_time),len(processed_freq)))
        for t in range(0,len(full_time)):
            min_time = 0
            for tg in range(0,len(gsm_time)):
	        if (float(full_time_sidereal[t])-float(gsm_time[tg]))>0:
                    min_time = tg 
                elif float(full_time_sidereal[t])==float(gsm_time[tg]):
                    min_time = tg       
            if float(full_time_sidereal[t])>float(gsm_time[-1]):
                min_time = -1
            for f in range(0,len(processed_freq)):
                min_freq = 0
                for fg in range(0,len(gsm_freq)):
                    if (float(processed_freq[f])-float(gsm_freq[fg]))>0:
                        min_freq = fg
                    elif float(processed_freq[f])==float(gsm_freq[fg]):
                        min_freq = fg
                        max_freq = fg
                if float(processed_freq[f])<=float(gsm_freq[0]):
                    min_freq = 0
                elif float(processed_freq[f])>=float(gsm_freq[-1]):
                    min_freq = -1
	        if full_data[t,f]>0.0:
                    if processed_freq[f]<60.:
                        full_gsm_data[t,f]=0.0
                    elif processed_freq[f]>90.:
                        full_gsm_data[t,f]=0.0 
                    else:
                        full_gsm_data[t,f] = gsm_data[min_time,min_freq]

#Full Data Plots (Including GSM)
if full:
    pylab.imshow(full_data,vmin=0,vmax=2e4,aspect=60./(full_time[-1]-full_time[0]),extent=(50,110,full_time[-1],full_time[0]))
    cbar = pylab.colorbar()
    cbar.set_label('Temperature (Kelvin)')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
    pylab.title('Full Dataset Waterfall Plot')
    pylab.savefig(result_dir+'full_waterfall_avgcal',dpi=300)
    pylab.clf()

    if gsm:
        pylab.imshow(full_gsm_data,vmin=0,vmax=2e4,aspect=60./(full_time[-1]-full_time[0]),extent=(50,110,full_time[-1],full_time[0]))
        cbar = pylab.colorbar() 
        cbar.set_label('Temperature (Kelvin)') 
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
        pylab.title('GSM Model Data for full dataset time/freq')
        pylab.savefig(result_dir+'full_waterfall_gsm_avgcal',dpi=300) 
        pylab.clf() 

        pylab.imshow(full_data[:,gsm_lim_ind[0]:gsm_lim_ind[-1]]-full_gsm_data[:,gsm_lim_ind[0]:gsm_lim_ind[1]],vmin=-2e3,vmax=2e3,aspect=30./(full_time[-1]-full_time[0]),extent=(60.,89.,full_time[-1],full_time[0]))
        cbar = pylab.colorbar()
        cbar.set_label('Temperature Difference (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
        pylab.title('Full Dataset - GSM Model Waterfall Plot')
        pylab.savefig(result_dir+'full_waterfall_data-gsm_avgcal',dpi=300)
        pylab.clf()

if eff:
    pylab.imshow(Eff_sm_mean,vmin=0,vmax=2e4,aspect=60./(time_delay[-1]*1e9),extent=(50,110,time_delay[-1]*1e9,0))
    cbar = pylab.colorbar()
    cbar.set_label('Temperature (Kelvin)')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Cal Time Delay (ns)')
    pylab.title('Efficiency Corrected Data')
    pylab.savefig(result_dir+'efficiency_waterfall_avgcal',dpi=300)
    pylab.clf()

#Processed Data Plots
pylab.scatter(processed_time,processed_data[:,120],c='b',edgecolor='b',s=3)
pylab.xlabel('Time (Hours since Midnight June 01 GMT)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Single Frequency (%0.1f MHz) Time Variation' %processed_freq[120])
pylab.grid()
pylab.xlim(processed_time[0],processed_time[-1])
pylab.ylim(0,7e3)
pylab.savefig(result_dir+'single_freq_time_var_avgcal',dpi=300)
pylab.clf()

pylab.scatter(sidereal_hour,processed_data[:,120],c='b',edgecolor='b',s=3)
pylab.xlabel('Time (Sidereal Hour)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Single Frequency (%0.1f MHz) Time Variation Folded by Day' %processed_freq[120])
pylab.grid()
pylab.xlim(0,24)
pylab.ylim(0,5e3)
pylab.savefig(result_dir+'single_freq_time_var_fold_avgcal',dpi=300)
pylab.clf()

pylab.scatter(processed_freq,mean_data,c='b',edgecolor='b',s=3,label='Mean Data')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Mean Data Spectrum')
pylab.grid()
pylab.xlim(50,110)
pylab.ylim(0,2e4)
pylab.savefig(result_dir+'mean_freq_spectrum_avgcal',dpi=300)

if fitting:
    pylab.plot(processed_freq,fitfunc(processed_freq,p[0],p[1],p[2]),c='g',label='Power Law fit (50-110 MHz)')
    pylab.plot(processed_freq,fitfunc(processed_freq,p1[0],p1[1],p1[2]),c='r',label='Power Law fit (%0.1f-%0.1f MHz)' %(data_lim[0],data_lim[1]))
    pylab.plot(processed_freq,fitfunc2(processed_freq,q[0],q[1]),c='c',label='-2.5 Power Law fit (50-110 MHz)')
    pylab.plot(processed_freq,fitfunc2(processed_freq,q1[0],q1[1]),c='m',label='-2.5 Power Law fit (%0.1f-%0.1f MHz)' %(data_lim[0],data_lim[1]))
    pylab.legend()
    pylab.title('Mean Data Spectrum with Power Law Fits')
    pylab.savefig(result_dir+'mean_freq_spectrum_avgcal_fit1_june1-5',dpi=300) 
    pylab.clf()

    pylab.scatter(processed_freq,mean_data-fitfunc(processed_freq,p1[0],p1[1],p1[2]),c='b',edgecolor='b',s=3,label='%0.1f-%0.1f MHz Fit (-%0.1f power)' %(data_lim[0],data_lim[1],p1[1]))
    pylab.scatter(processed_freq,mean_data-fitfunc2(processed_freq,q1[0],q1[1]),c='g',edgecolor='g',s=3,label='%0.1f-%0.1f MHz Fit (-2.5 power)' %(data_lim[0],data_lim[1]))
    pylab.xlabel('Frequency (MHz')
    pylab.ylabel('Temperature (Kelvin)')
    pylab.xlim(50,100)
    pylab.ylim(100,-100)
    pylab.grid()
    pylab.legend()
    pylab.savefig(result_dir+'mean_freq_spectrum_fit1_resid_june1-5',dpi=300)
    pylab.clf()

pylab.clf()  

pylab.scatter(sidereal_hour,processed_data[:,100],c='b',edgecolor='b',s=3)
pylab.scatter(sidereal_hour,processed_data[:,125],c='g',edgecolor='g',s=3)
pylab.scatter(sidereal_hour,processed_data[:,150],c='r',edgecolor='r',s=3)
pylab.scatter(sidereal_hour,processed_data[:,175],c='c',edgecolor='c',s=3)
pylab.scatter(sidereal_hour,processed_data[:,200],c='m',edgecolor='m',s=3)
pylab.xlabel('Time (Sidereal Hour)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Multiple Frequency Time Variation Folded by Day')
pylab.grid()
pylab.xlim(0,24)
pylab.ylim(0,1e4)
pylab.legend(('%0.1f MHz' %processed_freq[100],'%0.1f MHz' %processed_freq[125],'%0.1f MHz' %processed_freq[150],'%0.1f MHz' %processed_freq[175],'%0.1f MHz' %processed_freq[200]))
pylab.savefig(result_dir+'multi_freq_time_var_fold_avgcal',dpi=300)
pylab.clf()



