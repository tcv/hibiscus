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
gsm = True #If want to do GSM model comparison
svd = False #If want to do SVD calculations
full = True #If want to expand waterfall to fill in gaps
stack = False #If want to stack data by sidereal time
eff = False #If want to play with efficiency and time delay to uneff data


#Read in processed data
#data_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/'
data_dir = '/home/tcv/lustre/processed_data_subavg/'
result_dir = '/home/tcv/lustre/data_plots_subavg/'
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_guadalupe.dat')
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

data = os.listdir(data_dir)

processed_data = []
processed_time = []
processed_freq = []
processed_volt = []
#dates = ['01','03','04','05']
dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
date_files = []
day_ends = []
day_ends.append(0)

if gsm:
    gsm_freq = arange(50,111,1)
    gsm_time = arange(0,24,24./288)
    gsm_data = zeros((len(gsm_time),len(gsm_freq)))
    for t in range(0,len(gsm_time)):
        for f in range(0,len(gsm_freq)):
            gsm_data[t,f] = gsm_raw_data[t*61+f,2]

for i in range(0,len(dates)):
    curr_date = 'June'+dates[i]
    print 'Date being processed is:',curr_date
    single_data = loadtxt(data_dir+curr_date+'_processed_data_subavgcal.txt')
    single_data = array(single_data)
    for k in range(0,len(single_data)):
	processed_data.append(single_data[k])
    single_time = loadtxt(data_dir+curr_date+'_processed_time_subavgcal.txt')
    for j in range(0,len(single_time)):
        processed_time.append(single_time[j])
    single_volt = loadtxt(data_dir+curr_date+'_processed_volt_subavgcal.txt')
    for l in range(0,len(single_volt)):
        processed_volt.append(single_volt[l])
    if len(processed_freq)<1:
        processed_freq = loadtxt(data_dir+curr_date+'_processed_freq_subavgcal.txt')
    print 'Current Size of Processed Data is:',shape(processed_data)    
    day_ends.append(len(processed_data))
	    
processed_data = array(processed_data)
processed_time = array(processed_time)
print 'Full Shape of Processed Data is:',shape(processed_data)
print 'Shape of Processed Time is:',shape(processed_time)
print 'Nan Data present',where(isnan(processed_data))[0]
print 'Inf Data Present',where(isinf(processed_data))[0]

#Limited Frequencies to Try
data_lim = [65.,80.]
#data_lim = [65.,85.]
##################################################################################
data_lim_ind = [where(processed_freq<data_lim[0])[0][-1],where(processed_freq>data_lim[1])[0][0]]

# SVD for raw data:
if svd:
    U,S,Vh = linalg.svd(processed_data[:,data_lim_ind[0]:data_lim_ind[1]])
    #print U.shape,S.shape,Vh.shape
    Sm = linalg.diagsvd(S,len(processed_data),
                        len(processed_data[0,data_lim_ind[0]:data_lim_ind[1]]))
    S0 = S
    S0[0] = 0. 
    Sm0 = linalg.diagsvd(S0,len(processed_data),
                         len(processed_data[0,data_lim_ind[0]:data_lim_ind[1]]))
    P1 = dot(U,dot(Sm0,Vh))
    S0[1] = 0.
    Sm1 = linalg.diagsvd(S0,len(processed_data),
                         len(processed_data[0,data_lim_ind[0]:data_lim_ind[1]]))
    P2 =  dot(U,dot(Sm1,Vh))
    S0[2] = 0.
    Sm2 = linalg.diagsvd(S0,len(processed_data),
                         len(processed_data[0,data_lim_ind[0]:data_lim_ind[1]]))
    P3 =  dot(U,dot(Sm2,Vh))

# Adding Masking
prelim_mask = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(prelim_mask)):
    spike_mask = fc.spike_flag(processed_data[i],100)
    for j in range(0,len(prelim_mask[0])):
	if spike_mask[j] == 1.0:
	    prelim_mask[i,j] = 1.0
	    processed_data[i,j] = 0.0
        if processed_data[i,j]==0.0:
            prelim_mask[i,j] = 1.0

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
    mean_data[index]=0.0

#IF need to add efficiency data to other date. 
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic)
#mean_data = fc.effcal(Eff_sm,processed_freq,mean_data)

nandata = where(isnan(mean_data))
nandata = array(nandata)
for i in range(0,len(nandata[0])):
    index=nandata[0,i]
    mean_data[index]=0.0

#lgdata = where(mean_data>2e4)
#lgdata = array(lgdata)
#for i in range(0,len(lgdata[0])):
#    index = lgdata[0,i]

#Adding efficiency component
if eff:
    time_delay = arange(0,1e-8,5e-10)
    Eff_sm_data = []
    R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
    for i in range(0,len(time_delay)):
        R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,time_delay[i])
        Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
        Eff_sm = itp.UnivariateSpline(F_ant,Effic)
        new_mean_data = fc.effcal(Eff_sm,processed_freq,mean_data)
        Eff_sm_data.append(new_mean_data)
    

day_means = []
for i in range(0,len(day_ends)-1):
    single_day = []
    for j in range(0,len(processed_data[0])):
        single_freq = ma.array(processed_data[day_ends[i]:day_ends[i+1],j],mask=prelim_mask[day_ends[i]:day_ends[i+1],j])
        single_compress = ma.compressed(single_freq)
        single_mean = ma.mean(single_compress)
        single_day.append(single_mean)
    nandata = where(isnan(single_day))
    nandata = array(nandata)
    for k in range(0,len(nandata[0])):
        index = nandata[0,k]
        single_day[index] = 0.0
    day_means.append(single_day)

print 'Daily Mean Shape is:', shape(day_means)
day_means = array(day_means)

#Base Level power law fit
fitfunc = lambda x,p0,p1,p2: p0*x**(-p1)+p2
#rrfunc = lambda p,x,y: fitfunc(p,x)-y
pinit = [1e7,2.5,1e2]
p,success = opt.curve_fit(fitfunc,processed_freq,mean_data,pinit[:],ones(len(mean_data)))
#p,success = opt.leastsq(errfunc,p0,args=(processed_freq,mean_data))
p1,success1 = opt.curve_fit(fitfunc,processed_freq[data_lim_ind[0]:data_lim_ind[1]],mean_data[data_lim_ind[0]:data_lim_ind[1]],pinit[:],ones(len(mean_data[data_lim_ind[0]:data_lim_ind[1]])),maxfev=100000)
#p1,success1 = opt.leastsq(errfunc,p0,args=(processed_freq[data_lim_ind[0]:data_lim_ind[1]],mean_data[data_lim_ind[0]:data_lim_ind[1]]))
fitfunc2 = lambda x,q0,q1: q0*x**(-2.5)+q1
#errfunc2 = lambda q,a,b: fitfunc2(q,a)-b
qinit = [1e7,1e7]
q,success2 = opt.curve_fit(fitfunc2,processed_freq,mean_data,qinit[:],ones(len(mean_data)))
q1,success3= opt.curve_fit(fitfunc2,processed_freq[data_lim_ind[0]:data_lim_ind[1]],mean_data[data_lim_ind[0]:data_lim_ind[1]],qinit[:],ones(len(mean_data[data_lim_ind[0]:data_lim_ind[1]])))
#q,success2 = opt.leastsq(errfunc2,q0,args=(processed_freq,mean_data))
#q1,success3 = opt.leastsq(errfunc2,q0,args=(processed_freq[data_lim_ind[0]:data_lim_ind[1]],mean_data[data_lim_ind[0]:data_lim_ind[1]]))

print 'Full Data Fit Params (variable power law) are:',p
print 'Full Data Fit Param is:',q
print 'Limited Freq Data Fit Params (variable power law) are:',p1
print 'Limited Freq Data Fit Param is:', q1

day_params = []
for i in range(0,len(day_ends)-1):
    ptest,ptest_cov = opt.curve_fit(fitfunc,processed_freq[data_lim_ind[0]:data_lim_ind[1]],day_means[i,data_lim_ind[0]:data_lim_ind[1]],pinit[:],ones(len(mean_data[data_lim_ind[0]:data_lim_ind[1]])),maxfev=100000)
    day_params.append(ptest)

day_params = array(day_params)

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

if stack:
    sidereal_sort = argsort(sidereal_hour)
    sidereal_sort_time = sort(sidereal_hour)
    sidereal_data = zeros((len(processed_data),len(processed_data[0])))
    for i in range(0,len(sidereal_sort)):
        sidereal_data[i] = processed_data[sidereal_sort[i]]

    stack_time = gsm_time
    stack_data = zeros((len(stack_time),len(processed_freq)))
    for i in range(0,len(stack_time)):
        sub_data = []
        for j in range(0,len(sidereal_sort_time)):
            if abs(stack_time[i]-sidereal_sort_time[j])<1./24.:
                sub_data.append(sidereal_data[j])
        stack_data[i] = ma.mean(sub_data,axis=0)
    nandata = where(isnan(stack_data))
    nandata = array(nandata)
    for i in range(0,len(nandata[0])):
        index=nandata[0,i]
        stack_data[index]=zeros(len(stack_data[0]))
    infdata = where(isinf(stack_data))
    infdata = array(infdata)
    for i in range(0,len(infdata[0])):
        index=infdata[0,i]
        stack_data[index]=zeros(len(stack_data[0]))

    if svd:
        Us,Ss,Vhs = linalg.svd(stack_data[:,data_lim_ind[0]:data_lim_ind[1]])
        Ssm = linalg.diagsvd(Ss,len(stack_data),
                             len(stack_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Ss0 = Ss
        Ss0[0] = 0.
        Ssm0 = linalg.diagsvd(Ss0,len(stack_data),
                              len(stack_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Ps1 = dot(Us,dot(Ssm0,Vhs))
        Ss0[1] = 0.
        Ssm1 = linalg.diagsvd(Ss0,len(stack_data),
                              len(stack_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Ps2 = dot(Us,dot(Ssm1,Vhs))
        Ss0[2] = 0.
        Ssm2 = linalg.diagsvd(Ss0,len(stack_data),
                              len(stack_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Ps3 = dot(Us,dot(Ssm2,Vhs))

if full:
    full_time = arange(processed_time[0],processed_time[-1],0.05)
    full_data = zeros((len(full_time),len(processed_freq)))
    for i in range(0,len(full_data)):
        for j in range(0,len(processed_time)):
            if abs(full_time[i]-processed_time[j])<=0.1:
                full_data[i]=processed_data[j]
    
    if svd:
        Uf,Sf,Vhf = linalg.svd(full_data[:,data_lim_ind[0]:data_lim_ind[1]])
        Sfm = linalg.diagsvd(Sf,len(full_data),
                             len(full_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Sf0 = Sf
        Sf0[0] = 0.
        Sfm0 = linalg.diagsvd(Sf0,len(full_data),
                              len(full_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Pf1 = dot(Uf,dot(Sfm0,Vhf))
        Sf0[1] = 0.
        Sfm1 = linalg.diagsvd(Sf0,len(full_data),
                              len(full_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Pf2 = dot(Uf,dot(Sfm1,Vhf))
        Sf0[2] = 0.
        Sfm2 = linalg.diagsvd(Sf0,len(full_data),
                              len(full_data[0,data_lim_ind[0]:data_lim_ind[1]]))
        Pf3 = dot(Uf,dot(Sfm2,Vhf))

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
                    if processed_freq[f]<50.:
                        full_gsm_data[t,f]=0.0
                    elif processed_freq[f]>110.:
                        full_gsm_data[t,f]=0.0 
                    else:
                        full_gsm_data[t,f] = gsm_data[min_time,min_freq]

#Full Data Plots (Including GSM)
if full:
    pylab.imshow(full_data,vmin=0,vmax=2e4,aspect=100./(full_time[-1]-full_time[0]),extent=(40,140,full_time[-1],full_time[0]))
    cbar = pylab.colorbar()
    cbar.set_label('Temperature (Kelvin)')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
    pylab.title('Full Dataset Waterfall Plot')
    pylab.savefig(result_dir+'full_waterfall_avgcal',dpi=300)
    pylab.clf()

    if gsm:
        pylab.imshow(full_gsm_data,vmin=0,vmax=2e4,aspect=100./(full_time[-1]-full_time[0]),extent=(40,140,full_time[-1],full_time[0]))
        cbar = pylab.colorbar() 
        cbar.set_label('Temperature (Kelvin)') 
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
        pylab.title('GSM Model Data for full dataset time/freq')
        pylab.savefig(result_dir+'full_waterfall_gsm_avgcal',dpi=300) 
        pylab.clf() 

        pylab.imshow(full_data[:,40:287]-full_gsm_data[:,40:287],vmin=-2e3,vmax=2e3,aspect=60./(full_time[-1]-full_time[0]),extent=(50.,110.,full_time[-1],full_time[0]))
        cbar = pylab.colorbar()
        cbar.set_label('Temperature Difference (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
        pylab.title('Full Dataset - GSM Model Waterfall Plot')
        pylab.savefig(result_dir+'full_waterfall_data-gsm_avgcal',dpi=300)
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
pylab.xlim(40,140)
pylab.ylim(0,2e4)
pylab.savefig(result_dir+'mean_freq_spectrum_avgcal',dpi=300)
pylab.plot(processed_freq,fitfunc(processed_freq,p[0],p[1],p[2]),c='g',label='Power Law fit (40-140 MHz)')
pylab.plot(processed_freq,fitfunc(processed_freq,p1[0],p1[1],p1[2]),c='r',label='Power Law fit (%0.1f-%0.1f MHz)' %(data_lim[0],data_lim[1]))
pylab.plot(processed_freq,fitfunc2(processed_freq,q[0],q[1]),c='c',label='-2.5 Power Law fit (40-140 MHz)')
pylab.plot(processed_freq,fitfunc2(processed_freq,q1[0],q1[1]),c='m',label='-2.5 Power Law fit (%0.1f-%0.1f MHz)' %(data_lim[0],data_lim[1]))
pylab.legend()
pylab.title('Mean Data Spectrum with Power Law Fits')
pylab.savefig(result_dir+'mean_freq_spectrum_avgcal_fit',dpi=300) 
pylab.clf()

pylab.scatter(processed_freq,mean_data-fitfunc(processed_freq,p1[0],p1[1],p1[2]),c='b',edgecolor='b',s=3,label='60-80 MHz Fit (-%0.1f power)' %p1[1])
pylab.scatter(processed_freq,mean_data-fitfunc2(processed_freq,q1[0],q1[1]),c='g',edgecolor='g',s=3,label='60-80 MHz Fit (-2.5 power)')
pylab.xlabel('Frequency (MHz')
pylab.ylabel('Temperature (Kelvin)')
pylab.xlim(50,100)
pylab.ylim(100,-100)
pylab.grid()
pylab.legend()
pylab.savefig(result_dir+'mean_freq_spectrum_fit_resid',dpi=300)
pylab.clf()

pylab.scatter(processed_freq,day_means[0],c='b',edgecolor='b',s=1,label='June01')
pylab.scatter(processed_freq,day_means[2],c='g',edgecolor='g',s=1,label='June03')
pylab.scatter(processed_freq,day_means[3],c='c',edgecolor='c',s=1,label='June04')
pylab.scatter(processed_freq,day_means[4],c='r',edgecolor='r',s=1,label='June05')
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[0,0],day_params[0,1],day_params[0,2]),c='b',label='June01: -%0.2f Fit' %day_params[0,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[2,0],day_params[2,1],day_params[2,2]),c='g',label='June03: -%0.2f Fit' %day_params[2,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[3,0],day_params[3,1],day_params[3,2]),c='c',label='June04: -%0.2f Fit' %day_params[3,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[4,0],day_params[4,1],day_params[4,2]),c='r',label='June05: -%0.2f Fit' %day_params[4,1])
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Mean Daily Data Spectra and Fits')
pylab.grid()
pylab.xlim(40,140)
pylab.ylim(0,2e4)
pylab.legend()
pylab.savefig(result_dir+'mean_freq_spectrum_day_var_avgcal_fit_pt1',dpi=300)
pylab.clf()

pylab.scatter(processed_freq,day_means[5],c='y',edgecolor='y',s=1,label='June06')
pylab.scatter(processed_freq,day_means[7],c='m',edgecolor='m',s=1,label='June08')
pylab.scatter(processed_freq,day_means[8],c='k',edgecolor='k',s=1,label='June09')
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[5,0],day_params[5,1],day_params[5,2]),c='y',label='June06: -%0.2f Fit' %day_params[5,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[7,0],day_params[7,1],day_params[7,2]),c='m',label='June08: -%0.2f Fit' %day_params[7,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[8,0],day_params[8,1],day_params[8,2]),c='k',label='June09: -%0.2f Fit' %day_params[8,1])
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Temperature (Kelvin)') 
pylab.title('Mean Daily Data Spectra and Fits') 
pylab.grid() 
pylab.xlim(40,140) 
pylab.ylim(0,2e4) 
pylab.legend()
pylab.savefig(result_dir+'mean_freq_spectrum_day_var_avgcal_fit_pt2',dpi=300)
pylab.clf()

pylab.scatter(processed_freq,day_means[10],c='b',edgecolor='b',s=1,label='June11')
pylab.scatter(processed_freq,day_means[11],c='g',edgecolor='g',s=1,label='June12')
pylab.scatter(processed_freq,day_means[12],c='c',edgecolor='c',s=1,label='June13')
pylab.scatter(processed_freq,day_means[13],c='r',edgecolor='r',s=1,label='June14') 
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[10,0],day_params[10,1],day_params[10,2]),c='b',label='June11: -%0.2f Fit' %day_params[10,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[11,0],day_params[11,1],day_params[11,2]),c='g',label='June12: -%0.2f Fit' %day_params[11,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[12,0],day_params[12,1],day_params[12,2]),c='c',label='June13: -%0.2f Fit' %day_params[12,1])
pylab.plot(processed_freq,fitfunc(processed_freq,day_params[13,0],day_params[13,1],day_params[13,2]),c='r',label='June14: -%0.2f Fit' %day_params[13,1])
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Temperature (Kelvin)')  
pylab.title('Mean Daily Data Spectra and Fits')
pylab.grid()  
pylab.xlim(40,140)  
pylab.ylim(0,2e4)  
pylab.legend()
pylab.savefig(result_dir+'mean_freq_spectrum_day_var_avgcal_fit_pt3',dpi=300) 
pylab.clf() 

pylab.scatter(processed_freq,day_means[0]-fitfunc(processed_freq,day_params[0,0],day_params[0,1],day_params[0,2]),c='b',edgecolor='b',s=1,label='June01: -%0.2f Fit' %day_params[0,1])
pylab.scatter(processed_freq,day_means[2]-fitfunc(processed_freq,day_params[2,0],day_params[2,1],day_params[2,2]),c='g',edgecolor='g',s=1,label='June03: -%0.2f Fit' %day_params[2,1])
pylab.scatter(processed_freq,day_means[3]-fitfunc(processed_freq,day_params[3,0],day_params[3,1],day_params[3,2]),c='c',edgecolor='c',s=1,label='June04: -%0.2f Fit' %day_params[3,1])
pylab.scatter(processed_freq,day_means[4]-fitfunc(processed_freq,day_params[4,0],day_params[4,1],day_params[4,2]),c='r',edgecolor='r',s=1,label='June05: -%0.2f Fit' %day_params[4,1])
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Temperature (Kelvin)') 
pylab.title('Daily Mean Power Law Fit Residuals (Fit to %0.1f-%0.1fMHz)'%(data_lim[0],data_lim[1]))
pylab.grid()   
pylab.xlim(50,100)   
pylab.ylim(-100,100)   
pylab.legend()  
pylab.savefig(result_dir+'mean_freq_spectrum_day_var_fit_resid',dpi=300)
pylab.clf()  

pylab.scatter(processed_freq,day_means[5]-fitfunc(processed_freq,day_params[5,0],day_params[5,1],day_params[5,2]),c='y',edgecolor='y',s=1,label='June06: -%0.2f Fit' %day_params[5,1])
pylab.scatter(processed_freq,day_means[7]-fitfunc(processed_freq,day_params[7,0],day_params[7,1],day_params[7,2]),c='m',edgecolor='m',s=1,label='June08: -%0.2f Fit' %day_params[7,1])
pylab.scatter(processed_freq,day_means[8]-fitfunc(processed_freq,day_params[8,0],day_params[8,1],day_params[8,2]),c='k',edgecolor='k',s=1,label='June09: -%0.2f Fit' %day_params[8,1])
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Daily Mean Power Law Fit Residuals (Fit to %0.1f-%0.1fMHz)'%(data_lim[0],data_lim[1]))
pylab.grid()  
pylab.xlim(50,100)  
pylab.ylim(-100,100)  
pylab.legend() 
pylab.savefig(result_dir+'mean_freq_spectrum_day_var_fit_resid_pt2',dpi=300)
pylab.clf() 
 
pylab.scatter(processed_freq,day_means[10]-fitfunc(processed_freq,day_params[10,0],day_params[10,1],day_params[10,2]),c='b',edgecolor='b',s=1,label='June11: -%0.2f Fit' %day_params[10,1])
pylab.scatter(processed_freq,day_means[11]-fitfunc(processed_freq,day_params[11,0],day_params[11,1],day_params[11,2]),c='g',edgecolor='g',s=1,label='June12: -%0.2f Fit' %day_params[11,1])
pylab.scatter(processed_freq,day_means[12]-fitfunc(processed_freq,day_params[12,0],day_params[12,1],day_params[12,2]),c='c',edgecolor='c',s=1,label='June13: -%0.2f Fit')
pylab.scatter(processed_freq,day_means[13]-fitfunc(processed_freq,day_params[13,0],day_params[13,1],day_params[13,2]),c='r',edgecolor='r',s=1,label='June14: -%0.2f Fit')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Daily Mean Power Law Fit Residuals (Fit to %0.1f-%0.1fMHz)' %(data_lim[0],data_lim[1]))
pylab.grid()   
pylab.xlim(50,100)   
pylab.ylim(-100,100)   
pylab.legend() 
pylab.savefig(result_dir+'mean_freq_spectrum_day_var_fit_resid_pt3',dpi=300)
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


### SVD Plots
if svd:
    pylab.imshow(P1,vmin=-1e2,vmax=3e2,aspect=(data_lim[1]-data_lim[0])/(processed_time[-1]-processed_time[0]),extent=(data_lim[0],data_lim[1],processed_time[-1],processed_time[0]))
    cbar= pylab.colorbar()
    cbar.set_label('Temperature (Kelvin)')
    pylab.xlabel('Frequency (MHz)') 
    pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
    pylab.title('Dataset Waterfall Plot with 1 SVD Mode Removed')
    pylab.savefig(result_dir+'Unfill_waterfall_1SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
    pylab.clf() 

    pylab.imshow(P2,vmin=-1e2,vmax=3e2,aspect=(data_lim[1]-data_lim[0])/(processed_time[-1]-processed_time[0]),extent=(data_lim[0],data_lim[1],processed_time[-1],processed_time[0]))
    cbar= pylab.colorbar() 
    cbar.set_label('Temperature (Kelvin)') 
    pylab.xlabel('Frequency (MHz)')  
    pylab.ylabel('Time (Hours since Midnight June 01 GMT)')  
    pylab.title('Dataset Waterfall Plot with 2 SVD Modes Removed') 
    pylab.savefig(result_dir+'Unfill_waterfall_2SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300)
    pylab.clf()  

    pylab.imshow(P3,vmin=-1e2,vmax=3e2,aspect=(data_lim[1]-data_lim[0])/(processed_time[-1]-processed_time[0]),extent=(data_lim[0],data_lim[1],processed_time[-1],processed_time[0]))
    cbar= pylab.colorbar()  
    cbar.set_label('Temperature (Kelvin)')  
    pylab.xlabel('Frequency (MHz)')   
    pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
    pylab.title('Dataset Waterfall Plot with 3 SVD Modes Removed')  
    pylab.savefig(result_dir+'Unfill_waterfall_3SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300)
    pylab.clf()   

    pylab.imshow(processed_data[:,data_lim_ind[0]:data_lim_ind[1]]-dot(U,dot(Sm,Vh)),vmin=-1,vmax=1,aspect=(data_lim[1]-data_lim[0])/(processed_time[-1]-processed_time[0]),extent=(data_lim[0],data_lim[1],processed_time[-1],processed_time[0]))
    cbar= pylab.colorbar()  
    cbar.set_label('Temperature (Kelvin)')  
    pylab.xlabel('Frequency (MHz)')   
    pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
    pylab.title('Dataset - Full SVD (should be zeros)')  
    pylab.savefig(result_dir+'Unfill_waterfall_data-AllSVD_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300)
    pylab.clf()   

    if full:
        pylab.imshow(Pf1,vmin=-1e2,vmax=3e2,aspect=(data_lim[1]-data_lim[0])/(full_time[-1]-full_time[0]),extent=(data_lim[0],data_lim[1],full_time[-1],full_time[0]))
        cbar= pylab.colorbar() 
        cbar.set_label('Temperature (Kelvin)') 
        pylab.xlabel('Frequency (MHz)')  
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
        pylab.title('Full Dataset Waterfall with 1 SVD Mode Removed') 
        pylab.savefig(result_dir+'full_waterfall_1SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300)
        pylab.clf()  
 
        pylab.imshow(Pf2,vmin=-1e2,vmax=3e2,aspect=(data_lim[1]-data_lim[0])/(full_time[-1]-full_time[0]),extent=(data_lim[0],data_lim[1],full_time[-1],full_time[0]))
        cbar= pylab.colorbar()  
        cbar.set_label('Temperature (Kelvin)')  
        pylab.xlabel('Frequency (MHz)')   
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
        pylab.title('Full Dataset Waterfall with 2 SVD Modes Removed')  
        pylab.savefig(result_dir+'full_waterfall_2SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300)
        pylab.clf()   

        pylab.imshow(Pf3,vmin=-1e2,vmax=3e2,aspect=(data_lim[1]-data_lim[0])/(full_time[-1]-full_time[0]),extent=(data_lim[0],data_lim[1],full_time[-1],full_time[0]))
        cbar= pylab.colorbar()   
        cbar.set_label('Temperature (Kelvin)')   
        pylab.xlabel('Frequency (MHz)')    
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
        pylab.title('Full Dataset Waterfall with 3 SVD Modes Removed')   
        pylab.savefig(result_dir+'full_waterfall_3SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()    
 
        pylab.imshow(full_data[:,data_lim_ind[0]:data_lim_ind[1]]-dot(Uf,dot(Sfm,Vhf)),vmin=-1,vmax=1,aspect=(data_lim[1]-data_lim[0])/(full_time[-1]-full_time[0]),extent=(data_lim[0],data_lim[1],full_time[-1],full_time[0]))
        cbar= pylab.colorbar()   
        cbar.set_label('Temperature (Kelvin)')
        pylab.xlabel('Frequency (MHz)')  
        pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
        pylab.title('Full Dataset -Full SVD (should be zeros)')
        pylab.savefig(result_dir+'full_waterfall_data-AllSVD_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()

    if stack:
        pylab.imshow(Ps1,vmin=-5e2,vmax=5e2,aspect=(data_lim[1]-data_lim[0])/24.,extent=(data_lim[0],data_lim[1],24,0))
        cbar= pylab.colorbar()
        cbar.set_label('Temperature (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Sidereal Hour)')
        pylab.title('Stack Dataset Waterfall with 1 SVD Mode Removed')
        pylab.savefig(result_dir+'full_waterfall_stack_1SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()
  
        pylab.imshow(Ps2,vmin=-50,vmax=50,aspect=(data_lim[1]-data_lim[0])/24.,extent=(data_lim[0],data_lim[1],24,0))
        cbar= pylab.colorbar()
        cbar.set_label('Temperature (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Sidereal Hour)')
        pylab.title('Stack Dataset Waterfall with 2 SVD Modes Removed')
        pylab.savefig(result_dir+'full_waterfall_stack_2SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()
 
        pylab.imshow(Ps3,vmin=-50,vmax=50,aspect=(data_lim[1]-data_lim[0])/24.,extent=(data_lim[0],data_lim[1],24,0))
        cbar= pylab.colorbar()
        cbar.set_label('Temperature (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Sidereal Hour)')
        pylab.title('Stack Dataset Waterfall with 3 SVD Modes Removed')
        pylab.savefig(result_dir+'full_waterfall_stack_3SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()
  
        pylab.imshow(stack_data[:,data_lim_ind[0]:data_lim_ind[1]]-dot(Us,dot(Ssm,Vhs)),vmin=-1,vmax=1,aspect=(data_lim[1]-data_lim[0])/24.,extent=(data_lim[0],data_lim[1],24,0))
        cbar= pylab.colorbar()
        cbar.set_label('Temperature (Kelvin)')
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Time (Sidereal Hour)')
        pylab.title('Stack Dataset - Full SVD (should be zero)')
        pylab.savefig(result_dir+'full_waterfall_stack_data-AllSVD_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()

if svd:
    svd1_mean_data = []
    svd2_mean_data = []
    svd3_mean_data = []
    for i in range(0,len(P1[0])):
#    single_freq = ma.array(processed_data[:,i],mask=prelim_mask[:,i])
        single_freq_1svd = ma.array(P1[:,i],mask=prelim_mask[:,i])
        single_freq_2svd = ma.array(P2[:,i],mask=prelim_mask[:,i]) 
        single_freq_3svd = ma.array(P3[:,i],mask=prelim_mask[:,i])  
        single_compress_1svd = ma.compressed(single_freq_1svd)
        single_compress_2svd = ma.compressed(single_freq_2svd)
        single_compress_3svd = ma.compressed(single_freq_3svd)
#    single_compress = ma.compressed(single_freq)
#    single_mean = ma.mean(single_compress)
        single_mean_1svd = ma.mean(single_compress_1svd)
        svd1_mean_data.append(single_mean_1svd)
        single_mean_2svd = ma.mean(single_compress_2svd)
        svd2_mean_data.append(single_mean_2svd)
        single_mean_3svd = ma.mean(single_compress_3svd)
        svd3_mean_data.append(single_mean_3svd)
    nandata1 = where(isnan(svd1_mean_data))
    nandata1 = array(nandata1)
    for i in range(0,len(nandata1[0])):
        index=nandata1[0,i]
        svd1_mean_data[index]=0.0
    nandata2 = where(isnan(svd2_mean_data))
    nandata2 = array(nandata2)
    for i in range(0,len(nandata2[0])):
        index=nandata2[0,i]
        svd2_mean_data[index]=0.0
    nandata3 = where(isnan(svd3_mean_data))
    nandata3 = array(nandata3)
    for i in range(0,len(nandata3[0])):
        index=nandata3[0,i]
        svd3_mean_data[index]=0.0

    pylab.scatter(processed_freq[data_lim_ind[0]:data_lim_ind[1]],svd1_mean_data,c='b',edgecolor='b',s=3)
    pylab.scatter(processed_freq[data_lim_ind[0]:data_lim_ind[1]],svd2_mean_data,c='g',edgecolor='g',s=3)
    pylab.scatter(processed_freq[data_lim_ind[0]:data_lim_ind[1]],svd3_mean_data,c='r',edgecolor='r',s=3)
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Temperature (Kelvin)')
    pylab.title('Mean Data Spectrum')
    pylab.grid()
    pylab.xlim(data_lim[0],data_lim[1])
    pylab.ylim(-50,50)
    pylab.legend(('1 SVD Mode Rm','2 SVD Mode Rm','3 SVD Mode Rm'))
    pylab.savefig(result_dir+'mean_freq_spectrum_SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
    pylab.clf()

if stack:
    pylab.imshow(stack_data,vmin=0,vmax=2e4,aspect=100./24.,extent=(40.,140.,24.,0.))
    cbar = pylab.colorbar()
    cbar.set_label('Temperature (Kelvin)')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Time (Sidereal Hour)')
    pylab.title('Full Dataset Waterfall Folded on Sidereal Hour')
    pylab.savefig(result_dir+'full_waterfall_folded_avgcal',dpi=300)
    pylab.clf()

    if svd:
        svd1_mean_data = []
        svd2_mean_data = []
        svd3_mean_data = []
        for i in range(0,len(Ps1[0])):
            single_mean_1svd = ma.mean(Ps1[:,i])
            svd1_mean_data.append(single_mean_1svd)
            single_mean_2svd = ma.mean(Ps2[:,i])
            svd2_mean_data.append(single_mean_2svd)
            single_mean_3svd = ma.mean(Ps3[:,i])
            svd3_mean_data.append(single_mean_3svd)
        nandata1 = where(isnan(svd1_mean_data))
        nandata1 = array(nandata1)
        for i in range(0,len(nandata1[0])):
            index=nandata1[0,i]
            svd1_mean_data[index]=0.0
        nandata2 = where(isnan(svd2_mean_data))
        nandata2 = array(nandata2)
        for i in range(0,len(nandata2[0])):
            index=nandata2[0,i]
            svd2_mean_data[index]=0.0
        nandata3 = where(isnan(svd3_mean_data))
        nandata3 = array(nandata3)
        for i in range(0,len(nandata3[0])):
            index=nandata3[0,i]
            svd3_mean_data[index]=0.0

        pylab.scatter(processed_freq[data_lim_ind[0]:data_lim_ind[1]],svd1_mean_data,c='b',edgecolor='b',s=3)
        pylab.scatter(processed_freq[data_lim_ind[0]:data_lim_ind[1]],svd2_mean_data,c='g',edgecolor='g',s=3)
        pylab.scatter(processed_freq[data_lim_ind[0]:data_lim_ind[1]],svd3_mean_data,c='r',edgecolor='r',s=3)
        pylab.xlabel('Frequency (MHz)')
        pylab.ylabel('Temperature (Kelvin)')
        pylab.title('Mean Data Spectrum')
        pylab.grid()
        pylab.xlim(data_lim[0],data_lim[1])
        pylab.ylim(-50,50)
        pylab.legend(('1 SVD Mode Rm','2 SVD Mode Rm','3 SVD Mode Rm'))
        pylab.savefig(result_dir+'mean_freq_spectrum_stack_SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.ylim(-5,5)
        pylab.savefig(result_dir+'mean_freq_spectrum_stack_SVDrm_avgcal_'+str(int(data_lim[0]))+'-'+str(int(data_lim[1]))+'MHz',dpi=300) 
        pylab.clf()


