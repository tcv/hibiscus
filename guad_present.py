from numpy import *
import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import ephem as eph

#data_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/'
data_dir = '/home/tcv/lustre/processed_data_take2/'
result_dir = '/home/tcv/lustre/data_plots_take2/'
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_guadalupe.dat')

data = os.listdir(data_dir)

processed_data = []
processed_time = []
processed_freq = []

dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
date_files = []

gsm_freq = arange(50,111,1)
gsm_time = arange(0,24,24./288)
gsm_data = zeros((len(gsm_time),len(gsm_freq)))
for t in range(0,len(gsm_time)):
    for f in range(0,len(gsm_freq)):
        gsm_data[t,f] = gsm_raw_data[t*61+f,2]

for i in range(0,len(dates)):
    curr_date = 'June'+dates[i]
    print 'Date being processed is:',curr_date
    single_data = loadtxt(data_dir+curr_date+'_processed_data_avgcal.txt')
    single_data = array(single_data)
    for k in range(0,len(single_data)):
	processed_data.append(single_data[k])
    single_time = loadtxt(data_dir+curr_date+'_processed_time_avgcal.txt')
    for j in range(0,len(single_time)):
        processed_time.append(single_time[j])
    if len(processed_freq)<1:
        processed_freq = loadtxt(data_dir+curr_date+'_processed_freq_avgcal.txt')
    print 'Current Size of Processed Data is:',shape(processed_data)    
	    
#    for file in data:
#	curr_date = 'June'+dates[i]
#	if file.split('-')[0]==curr_date:
#	    if file.split('.')[-1]=='txt':
#	        date_files.append(file)
#    print 'Number of Files in that day is:',len(date_files)/3
#    num_files = len(date_files)/3
#    for j in range(0,num_files):
#	full_date = curr_date+'-'+str(j)
#	if full_date != 'June09-5':
#      	    single_data = loadtxt(data_dir+full_date+'_processed_data_avgcal.txt')
# 	    single_data = array(single_data)
#	    if len(single_data)<550:
#	        for k in range(0,len(single_data)):
#		    processed_data.append(single_data[k])
#	    else:
#	        processed_data.append(single_data)
#	    single_time = loadtxt(data_dir+full_date+'_processed_time_avgcal.txt')
#	    processed_time.append(single_time)
#	    if len(processed_freq)<1:
#	        processed_freq = loadtxt(data_dir+full_date+'_processed_freq_avgcal.txt')
#	    print 'Current status of processed data:',shape(processed_data)
#    date_files = []


#print 'Full Shape of Processed Data is:',shape(processed_data)
processed_data = array(processed_data)
processed_time = array(processed_time)
print 'Full Shape of Processed Data is:',shape(processed_data)
print 'Shape of Processed Time is:',shape(processed_time)

#What plots are most useful?

prelim_mask = zeros((len(processed_data),len(processed_data[0])))

for i in range(0,len(prelim_mask)):
    spike_mask = fc.spike_flag(processed_data[i],100)
    for j in range(0,len(prelim_mask[0])):
	if spike_mask[j] == 1.0:
	    prelim_mask[i,j] = 1.0
	    processed_data[i,j] = 0.0
        if processed_data[i,j]==0.0:
            prelim_mask[i,j] = 1.0

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


initial = eph.date('2013/6/1')
ephem_dates = []
for i in range(0,len(processed_time)):
    ephem_dates.append(eph.date(initial+processed_time[i]/24.))
    

guad = eph.Observer()
guad.lon = '-118.3'
guad.lat = '28.8833'

ephem_sidereal = []
for i in range(0,len(processed_time)):
    guad.date = ephem_dates[i]
    test = guad.sidereal_time()
    ephem_sidereal.append(test)

sidereal_hour = []
for i in range(0,len(processed_time)):
    sidereal_hour.append(ephem_sidereal[i]*12./pi)

full_time = arange(processed_time[0],processed_time[-1],0.05)
full_data = zeros((len(full_time),len(processed_freq)))
for i in range(0,len(full_data)):
    for j in range(0,len(processed_time)):
        if abs(full_time[i]-processed_time[j])<=0.1:
            full_data[i]=processed_data[j]

full_time_sidereal = []
for t in range(0,len(full_time)):
    single_date = eph.date(initial+full_time[t]/24.)
    guad.date = single_date
    single_time = guad.sidereal_time()
    single_hour = single_time*12./pi
    full_time_sidereal.append(single_hour)

gsm_diff = zeros((len(full_time),len(processed_freq)))
for t in range(0,len(full_time)):
    min_time = 0
    max_time = 0
    for tg in range(0,len(gsm_time)):
	if full_time[t]<gsm_time[tg]:
            min_time = tg 
            max_time = tg+1
        elif full_time[t]==gsm_time[tg]:
            min_time = tg
            max_time = tg        
    if full_time[t]>gsm_time[-1]:
        min_time = -1
        max_time = -1
    for f in range(0,len(processed_freq)):
	min_freq = 0
        max_freq = 0
        for fg in range(0,len(gsm_freq)):
            if processed_freq[f]<gsm_freq[fg]:
                min_freq = fg
                max_freq = fg+1
            elif processed_freq[f]==gsm_freq[fg]:
                min_freq = fg
                max_freq = fg
        if processed_freq[f]<=gsm_freq[0]:
            min_freq = 0
            max_freq = 0
        elif processed_freq[f]>=gsm_freq[-1]:
            min_freq = -1
            max_freq = -1
	if full_data[t,f]>0.0:
#            if min_time == max_time:
#                if min_freq == max_freq:
#                    gsm_diff[t,f] = full_data[t,f]-gsm_data[min_time,min_freq]
#                else:
#                    (Ta,Tb) = polyfit(gsm_freq[min_freq:max_freq],gsm_data[min_time,min_freq:max_freq],1)
#                    gsm_diff[t,f] = full_data[t,f]-polyval([Ta,Tb],processed_freq[f])
#            else:
#                if min_freq == max_freq:
#                    (Ta,Tb) = polyfit(gsm_time[min_time:max_time],gsm_data[min_time:max_time,min_freq],1)
#                    gsm_diff[t,f] = full_data[t,f]-polyval([Ta,Tb],full_time[t]) 
#                else:
#		    (Tamin,Tbmin) = polyfit(gsm_freq[min_freq:max_freq],gsm_data[min_time,min_freq:max_freq],1)
#                    (Tamax,Tbmax) = polyfit(gsm_freq[min_freq:max_freq],gsm_data[max_time,min_freq:max_freq],1)
#                    tdata = [polyval([Tamin,Tbmin],processed_freq[f]),polyval([Tamax,Tbmax],processed_freq[f])]
#                    (Ta,Tb) = polyfit(gsm_time[min_time:max_time],tdata,1)
#                    gsm_diff[t,f] = full_data[t,f]-polyval([Ta,Tb],full_time[t])
            if processed_freq[f]<50.:
                gsm_diff[t,f]=0.0
            elif processed_freq[f]>110.:
                gsm_diff[t,f]=0.0 
            else:
                gsm_diff[t,f] = full_data[t,f]-gsm_data[min_time,min_freq]
          

pylab.imshow(full_data,vmin=0,vmax=2e4,aspect=100./(full_time[-1]-full_time[0]),extent=(40,140,full_time[-1],full_time[0]))
cbar = pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_avgcal',dpi=300)
pylab.clf()

pylab.imshow(gsm_diff,vmin=-5e2,vmax=5e3,aspect=100./(full_time[-1]-full_time[0]),extent=(40,140,full_time[-1],full_time[0]))
cbar = pylab.colorbar() 
cbar.set_label('Temperature (Kelvin)') 
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
pylab.title('Full Dataset Waterfall Plot - GSM Model Data') 
pylab.savefig(result_dir+'full_waterfall_subtract_avgcal',dpi=300) 
pylab.clf() 

pylab.scatter(processed_time,processed_data[:,120],c='b',edgecolor='b',s=3)
pylab.xlabel('Time (Hours since Midnight June 01 GMT)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Single Frequency (%0.1f MHz) Time Variation' %processed_freq[120])
pylab.grid()
pylab.xlim(processed_time[0],processed_time[-1])
pylab.ylim(0,5e3)
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

pylab.scatter(processed_freq,mean_data,c='b',edgecolor='b',s=3)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Mean Data Spectrum')
pylab.grid()
pylab.xlim(40,140)
pylab.ylim(0,2e4)
pylab.savefig(result_dir+'mean_freq_spectrum_avgcal',dpi=300)
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


