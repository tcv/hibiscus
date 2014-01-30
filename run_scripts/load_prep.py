"""
Module to calculate the 50 Ohm spectrum for a single day of data
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
import file_funcs as ff
import cal_funcs as cf


#Main directories for the input and output
indir = '/lustre/tcv/file_rebinned_data/'
outdir='/lustre/tcv/calibration_data/'
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
        if fname.split('_')[-1]=='50ohm.dat':
            filename = directory+fname
#load data file
            time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40.+width/2,130.,width)
            if len(freq)>len(sub_data):
                freq = freq[0:-2]
            load_data.append(sub_data)
            short_time.append(time)
        elif fname.split('_')[-1]=='mask.dat':
            if fname.split('_')[-2]=='50ohm':
                filename = directory+fname
#load mask file
                time,form,sub_mask,mask,freq,volt,temp = ff.loadsingle(filename)
                load_mask.append(sub_mask)
                load_mtime.append(time)        

load_data = array(load_data)
load_mask = array(load_mask)

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

load_array = ma.array(mean_load,mask=mean_mask)
load_comp = ma.compressed(load_array)
freq_array = ma.array(freq,mask=mean_mask)
freq_comp = ma.compressed(freq_array)

#Need to limit the frequencies for fit just to what we are considering.
fmin = where(freq_comp<=60.)[0][-1]
fmax = where(freq_comp<=90.)[0][-1]
(Fa,Fb) = polyfit(freq_comp[fmin:fmax],load_comp[fmin:fmax],1)

savetxt(outdir+'June_'+date_ind+'_avg_50ohm.txt',polyval([Fa,Fb],freq),delimiter=' ')

#Plotting Checks
pylab.imshow(sortload*10**9,aspect=90./len(sortind),extent=(40,130,len(sort_ind),0.0))
pylab.colorbar()
pylab.title('Variation of 50 Ohm over the day')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (50 Ohm index)')
pylab.savefig(outdir+'June_'+date_ind+'_50ohm_variation',dpi=300)
pylab.clf()

pylab.plot(freq_comp,load_comp*10**9,label='Mean Data')
pylab.plot(freq,polyval([Fa,Fb],freq)*10**9,label='Linear Fit')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Power (nW)')
pylab.grid()
pylab.tile('Mean 50 Ohm and Fit')
pylab.savefig(outdir+'June_'+date_ind+'_50ohm_mean',dpi=300)
pylab.clf()