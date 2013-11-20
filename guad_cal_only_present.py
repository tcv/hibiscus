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

data_dir = '/home/tcv/lustre/calibrated_data_test/'

processed_data = []
processed_time = []
processed_volt = []
processed_std = []
comp_day_data = []
comp_day_time = []
comp_day_volt = []
day_ends = []
day_ends.append(0)

dates = ['01','02','03','04','05','06','07','08','09','11','12','13','14']

testfile = loadtxt(data_dir+'June01_day_to_night/'+os.listdir(data_dir+'June01_day_to_night')[0])

freqs = arange(0,250.,250./len(testfile))
binds = []
einds = []
fmin = 45.
fint = 0.2
for i in range(0,255):
    bind = where(freqs<=(fmin+fint*i))[0][-1]
    eind = where(freqs>=(fmin+fint*(i+1)))[0][0]
#    print 'Sum in Bin:', eind-bind
    binds.append(bind)
    einds.append(eind)    

for date in dates:
    direct = 'June'+date+'_day_to_night'
    print direct
    directory = data_dir+direct+'/'
    dirlist = os.listdir(directory)
    for fname in dirlist:
        if fname.split('_')[-2]=='antenna':
            time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(directory+fname)
            nandata = where(isnan(sub_data))
            nandata = array(nandata)
            for i in range(0,len(nandata)):
                sub_data[nandata[i]] = 0.0
            mask = fc.flagging(sub_data,freqs,3.,131)
            smask = fc.spike_flag(sub_data,100)
            for i in range(0,len(smask)):
                if smask[i]==1.0:
                    mask[i] =1.0
            new_data = []
            new_std = []
            for i in range(0,len(binds)):
                single_array = ma.array(sub_data[binds[i]:einds[i]],mask=mask[binds[i]:einds[i]])
                single_comp = ma.compressed(single_array)
                single_mean = ma.mean(single_comp)
                single_std = ma.std(single_comp)
                new_data.append(single_mean)
                new_std.append(single_std)
            processed_data.append(new_data)
            processed_std.append(new_std)
            processed_time.append(time)
            processed_volt.append(volt)
            if date =='01':            
                comp_day_data.append(new_data)
                comp_day_time.append(time) 
                comp_day_volt.append(volt)
    day_ends.append(len(processed_data))


processed_data = array(processed_data)
processed_time = array(processed_time)
processed_volt = array(processed_volt)
processed_std = array(processed_std)
comp_day_data = array(comp_day_data)
comp_day_time = array(comp_day_time)
comp_day_volt = array(comp_day_volt)

mean_processed_data = ma.mean(processed_data,axis=0)
std_processed_data = ma.std(processed_data,axis=0)

mask = zeros((len(processed_data),len(processed_data[0])))
for i in range(1,len(processed_data)-1):
    for j in range(0,len(processed_data[0])):
        if processed_volt[i]>=12.:
            mask[i,j] = 1.0
        elif processed_volt[i]<=9.75:
            mask[i,j] = 1.0
#    for j in range(0,len(processed_data[0])):
        if processed_data[i,j]>(mean_processed_data[j]+2*std_processed_data[j]):
            mask[i,j] = 1.0
        elif processed_data[i,j]<(mean_processed_data[j]-2*std_processed_data[j]):
            mask[i,j] = 1.0
#    if processed_data[i]>(0.5*(processed_data[i-1]+processed_data[i+1])+500):
#        mask[i] = 1.0
        if isnan(processed_data[i,j])==True:
            mask[i,j] = 1.0
            processed_data[i,j] = 0.
        if isinf(processed_data[i,j])==True:
            mask[i,j] = 1.0
            processed_data[i] = 0.

print 'Percentage Masked Data:', 100*sum(mask,axis=0)/len(mask)
#print where(isnan(processed_data)), where(isinf(processed_data))

comp_mask = zeros((len(comp_day_data),len(comp_day_data[0])))
for i in range(0,len(comp_day_data)):
    for j in range(0,len(processed_data[0])):
        if comp_day_volt[i]>=12.:
            comp_mask[i,j] = 1.0
        elif comp_day_volt[i]<=9.75:
            comp_mask[i,j] = 1.0
        if comp_day_data[i,j]>(mean_processed_data[j]+2*std_processed_data[j]):
            comp_mask[i,j] = 1.0
        elif comp_day_data[i,j]<(mean_processed_data[j]-2*std_processed_data[j]):
            comp_mask[i,j] = 1.0
#    if comp_day_data[i]>(0.5*(comp_day_data[i-1]+comp_day_data[i+1])+500):
#        comp_mask[i] = 1.0
        if isnan(comp_day_data[i,j])==True: 
            comp_mask[i,j] = 1.0 
            comp_day_data[i,j] = 0. 
        if isinf(comp_day_data[i,j])==True: 
            comp_mask[i,j] = 1.0 
            comp_day_data[i,j] = 0. 


initial = eph.date('2013/6/1')
guad = eph.Observer()
guad.lon = '-118.3'
guad.lat = '28.8833'
sidereal_hour = []
comp_day_hour = []
for i in range(0,len(processed_time)):
    single_date = eph.date(initial+processed_time[i]/24.)
    guad.date = single_date
    single_time = guad.sidereal_time()
    sidereal_hour.append(single_time*12./pi)

for i in range(0,len(comp_day_data)):
    single_date = eph.date(initial+comp_day_time[i]/24.)
    guad.date = single_date
    single_time = guad.sidereal_time()
    comp_day_hour.append(single_time*12./pi)

comp_day_sort = argsort(comp_day_hour)
cdtime = sort(comp_day_hour)
cddata = zeros((len(comp_day_data),len(comp_day_data[0])))
cdmask = zeros((len(comp_day_data),len(comp_day_data[0])))
for i in range(0,len(comp_day_sort)):
    cddata[i] = comp_day_data[comp_day_sort[i]]
    cdmask[i] = comp_mask[comp_day_sort[i]]

refit_data = zeros((len(processed_data),len(processed_data[0])))
for i in range(1,len(day_ends)):
    diff_data = []
    for f in range(0,len(processed_data[0])):
#        diff_data = []
        for j in range(day_ends[i-1],day_ends[i]):
            ct = 24.
            ctind = 0
            for k in range(0,len(cddata)):
                if abs(sidereal_hour[j]-cdtime[k])<ct:
                    if cdmask[k,f]==0.:
                        ct = abs(sidereal_hour[j]-cdtime[k])
                        ctind = k 
#        for f in range(0,15):
            if mask[j,f]==0.:
#            if cdmask[ctind]==0.:
                diff_data.append(processed_data[j,f]-cddata[ctind,f])
    avg_diff = ma.mean(diff_data)
    print 'June ', dates[i-1], ' offset:',avg_diff
    for f in range(0,len(processed_data[0])):
        for j in range(day_ends[i-1],day_ends[i]):
            refit_data[j,f] = processed_data[j,f]-avg_diff 
 

sidereal_sort = argsort(sidereal_hour) 
stime = sort(sidereal_hour) 
sdata = zeros((len(processed_data),len(processed_data[0])))
smask = zeros((len(processed_data),len(processed_data[0])))
rfsdata = zeros((len(processed_data),len(processed_data[0])))
stddata = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(sidereal_sort)): 
    sdata[i] = processed_data[sidereal_sort[i]] 
    smask[i] = mask[sidereal_sort[i]] 
    rfsdata[i] = refit_data[sidereal_sort[i]]    
    stddata[i] = processed_std[sidereal_sort[i]]

for i in range(1,len(smask)-1):
    for f in range(0,len(processed_data[0])):
        if sdata[i,f]>=(0.5*(sdata[i-1,f]+sdata[i+1,f])+500):
            smask[i,f] = 1.0
        elif sdata[i,f]<=(0.5*(sdata[i-1,f]+sdata[i+1,f])-500):
            smask[i,f] = 1.0

print 'New Mask percentages:',100.*sum(smask,axis=0)/len(smask)

stack_time = arange(0,24,0.01)
stack_data = zeros((len(stack_time),len(processed_data[0])))
rfstack_data = zeros((len(stack_time),len(processed_data[0])))
stack_std = zeros((len(stack_time),len(processed_data[0])))
rfstack_std = zeros((len(stack_time),len(processed_data[0])))
max_ind = where(stime>=stack_time[1])[0][0]
for f in range(0,len(processed_data[0])):
    new_data = ma.array(sdata[0:max_ind,f],mask=smask[0:max_ind,f])
    new_std = ma.array(stddata[0:max_ind,f],mask=smask[0:max_ind,f])
    new_compress = ma.compressed(new_data)
    new_cstd = ma.compressed(new_std)
    cw = 1/new_cstd**2
    stack_data[0,f] = ma.average(new_compress,weights=cw)
    stack_std[0,f] = sqrt(sum(cw*(new_compress-stack_data[0,f]*ones(len(new_compress)))**2)/sum(cw))
    rf_new_data = ma.array(rfsdata[0:max_ind,f],mask=smask[0:max_ind,f])
    rf_new_compress = ma.compressed(rf_new_data)
    rfstack_data[0,f] = ma.average(rf_new_compress,weights=cw)
    rfstack_std[0,f] = sqrt(sum(cw*(rf_new_compress-rfstack_data[0,f]*ones(len(rf_new_compress)))**2)/sum(cw))
    for i in range(1,len(stack_time)-1):
        min_ind = where(stime<=stack_time[i])[0][-1]
        max_ind = where(stime>=stack_time[i+1])[0][0]
        new_data = ma.array(sdata[min_ind:max_ind,f],mask=smask[min_ind:max_ind,f])
        new_std = ma.array(stddata[min_ind:max_ind,f],mask=smask[min_ind:max_ind,f])
        new_compress = ma.compressed(new_data)
        new_cstd = ma.compressed(new_std)
        cw = 1/new_cstd**2
        stack_data[i,f] = ma.average(new_compress,weights=cw)
        stack_std[i,f] = sqrt(sum(cw*(new_compress-stack_data[i,f]*ones(len(new_compress)))**2)/sum(cw))
        rf_new_data = ma.array(rfsdata[min_ind:max_ind,f],mask=smask[min_ind:max_ind,f])
        rf_new_compress = ma.compressed(rf_new_data)
        rfstack_data[i,f] = ma.average(rf_new_compress,weights=cw)
        rfstack_std[i,f] = sqrt(sum(cw*(rf_new_compress-rfstack_data[i,f]*ones(len(rf_new_compress)))**2)/sum(cw))

    min_ind = where(stime<=stack_time[-1])[0][-1]
    new_data = ma.array(sdata[min_ind:-1,f],mask=smask[min_ind:-1,f])
    new_std = ma.array(stddata[min_ind:-1,f],mask=smask[min_ind:-1,f])
    new_compress = ma.compressed(new_data)
    new_cstd = ma.compressed(new_std)
    cw = 1/new_cstd**2
    stack_data[-1,f] = ma.average(new_compress,weights=cw)
    stack_std[-1,f] = sqrt(sum(cw*(new_compress-stack_data[-1,f]*ones(len(new_compress)))**2)/sum(cw)) 
    rf_new_data = ma.array(rfsdata[min_ind:-1,f],mask=smask[min_ind:-1,f])
    rf_new_compress = ma.compressed(rf_new_data) 
    rfstack_data[-1,f] = ma.average(rf_new_compress,weights=cw) 
    rfstack_std[-1,f] = sqrt(sum(cw*(rf_new_compress-rfstack_data[-1,f]*ones(len(rf_new_compress)))**2)/sum(cw)) 

savetxt('/home/tcv/stack_data.txt',stack_data,delimiter=' ')
savetxt('/home/tcv/stack_std.txt',stack_std,delimiter=' ')
savetxt('/home/tcv/rfstack_data.txt',rfstack_data,delimiter=' ')
savetxt('/home/tcv/rfstack_std.txt',rfstack_std,delimiter=' ')
savetxt('/home/tcv/stack_time.txt',stack_time,delimiter=' ')


#pylab.scatter(stack_time,stack_data,c='b',edgecolor='b',s=3)
#pylab.ylim(0,5000.) 
#pylab.xlim(0,24) 
#pylab.grid() 
#pylab.xlabel('Local Sidereal Time (Hours)') 
#pylab.ylabel('Temperature (Kelvin)') 
#pylab.legend(('70 MHz Sky Signal')) 
#pylab.savefig('/home/tcv/diurnal_stack_70MHz',dpi=300)
#pylab.clf()

#pylab.errorbar(stack_time,stack_data,stack_std)
#pylab.ylim(0,5000.)
#pylab.xlim(0,24)
#pylab.grid()
#pylab.xlabel('Local Sidereal Time (Hours)')
#pylab.ylabel('Temperature (Kelvin)')
#pylab.legend(('70 MHz Sky Signal'))
#pylab.savefig('/home/tcv/diurnal_stack_werror_70MHz',dpi=300)
#pylab.clf()

#pylab.scatter(stack_time,rfstack_data,c='b',edgecolor='b',s=3)
#pylab.ylim(0,5000.)  
#pylab.xlim(0,24)  
#pylab.grid()  
#pylab.xlabel('Local Sidereal Time (Hours)')
#pylab.ylabel('Temperature (Kelvin)')
#pylab.legend(('70 MHz Sky Signal')) 
#pylab.savefig('/home/tcv/diurnal_stack_daycorr_70MHz',dpi=300)
#pylab.clf() 
 
#pylab.errorbar(stack_time,rfstack_data,rfstack_std)
#pylab.ylim(0,5000.) 
#pylab.xlim(0,24) 
#pylab.grid() 
#pylab.xlabel('Local Sidereal Time (Hours)')
#pylab.ylabel('Temperature (Kelvin)')
#pylab.legend(('70 MHz Sky Signal'))
#pylab.savefig('/home/tcv/diurnal_stack_daycorr_werror_70MHz',dpi=300)
#pylab.clf() 
    
