"""
Module to calculate the noise spectrum for a single day of data
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
import os
import skrf as rf
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('/home/tcv/hibiscus'))
import file_funcs as ff
import cal_funcs as cf


#Main directories for the input and output
indir = '/lustre/tcv/freq_rebinned_data/'
outdir='/lustre/tcv/calibration_data/'
#indir = '/lustre/tcv/truncated_data/'
#outdir = '/lustre/tcv/rfi_check_data/'
directories = os.listdir(indir)

#Setting a single day for parallel computation
date_ind = sys.argv[1]

load_data = []
load_mask = []
load_time = []
load_mtime = []
if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being rebinned is:',direct
    directory = indir+direct+'/'
    new_directory = outdir+direct+'/'
    dirlist = os.listdir(directory)
#Iterate for each file in the directory
    for fname in dirlist:
        if fname.split('_')[-1]=='noise.dat':
            filename = directory+fname
#load data file
            time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40.+width/2,130.,width)
            if len(freq)>len(sub_data):
                freq = freq[0:-2]
            load_data.append(sub_data)
            load_time.append(time)
        elif fname.split('_')[-1]=='noise.dat':
            if fname.split('_')[-2]=='noise':
                filename = directory+fname
#load mask file
                time,form,sub_mask,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 90.0/len(sub_mask)
                freq = arange(40.+width/2,130.,width)
                if len(freq)>len(sub_mask):
                    freq = freq[0:-2]
                load_mask.append(sub_mask)
                load_mtime.append(time)        

load_data = array(load_data)
load_mask = array(load_mask)

if len(load_mask)==0:
    load_mask = zeros((len(load_data),len(load_data[0])))
    load_mtime = load_time

sortind = argsort(load_time)
sorttime = zeros(len(load_data))
sortload = zeros((len(load_data),len(load_data[0])))
sortindm = argsort(load_mtime)
sortmask = zeros((len(load_data),len(load_data[0])))
for i in range(0,len(sortind)):
    sorttime[i] = load_time[sortind[i]]
    sortload[i] = load_data[sortind[i]]
    sortmask[i] = load_mask[sortindm[i]]

mean_load, mean_mask = cf.time_mean(sortload,sortmask)
mean_mask = ff.spike_flag(mean_load,mean_mask,freq,2.)
print freq[where(mean_mask==1.0)[0]]

load_array = ma.array(mean_load,mask=mean_mask)
load_comp = ma.compressed(load_array)
freq_array = ma.array(freq,mask=mean_mask)
freq_comp = ma.compressed(freq_array)

#Need to limit the frequencies for fit just to what we are considering.
fmin = where(freq_comp<=50.)[0][-1]
fmax = where(freq_comp<=100.)[0][-1]
(Fa,Fb,Fc) = polyfit(freq_comp[fmin:fmax],load_comp[fmin:fmax],2)

savetxt(outdir+'June_'+date_ind+'_avg_noise.txt',mean_load,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_fit_noise.txt',polyval([Fa,Fb,Fc],freq),delimiter=' ')

#Plotting Checks
pylab.imshow(sortload*10**9,vmin=15,vmax=30,aspect=90./len(sortind),extent=(40,130,len(sortind),0.0))
pylab.colorbar()
#pylab.title('Variation of 100 Ohm over the day')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (100 Ohm index)')
pylab.savefig(outdir+'June_'+date_ind+'_noise_variation',dpi=300)
pylab.clf()

pylab.scatter(freq,mean_load*10**9,c='b',edgecolor='b',label='Unfiltered Mean Data')
pylab.plot(freq_comp,load_comp*10**9,label='Mean Data')
pylab.plot(freq,polyval([Fa,Fb,Fc],freq)*10**9,label='Linear Fit')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Power (nW)')
pylab.ylim(0,30)
pylab.xlim(40,130)
pylab.legend()
pylab.grid()
#pylab.title('Mean Noise and Fit')
pylab.savefig(outdir+'June_'+date_ind+'_noise_mean',dpi=300)
pylab.clf()
