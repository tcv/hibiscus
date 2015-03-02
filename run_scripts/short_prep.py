"""
Module to calculate the short spectrum for a single day of data
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
#indir = '/lustre/tcv/freq_rebinned_data/'
#outdir='/lustre/tcv/calibration_data/'
#directories = os.listdir(indir)
indir = '/lustre/tcv/truncated_data/'
outdir = '/lustre/tcv/rfi_check_data/'
directories = os.listdir(indir)

#Setting a single day for parallel computation
date_ind = sys.argv[1]

short_data = []
short_mask = []
short_time = []
short_mtime = []
if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being rebinned is:',direct
    directory = indir+direct+'/'
    new_directory = outdir+direct+'/'
    dirlist = os.listdir(directory)
#Iterate for each file in the directory
    for fname in dirlist:
        if fname.split('_')[-1]=='short.dat':
            filename = directory+fname
#load data file
            time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40.+width/2,130.,width)
            if len(freq)>len(sub_data):
                freq = freq[0:-2]
            short_data.append(sub_data)
            short_time.append(time)
        elif fname.split('_')[-1]=='mask.dat':
            if fname.split('_')[-2]=='short':
                filename = directory+fname
#load mask file
                time,form,sub_mask,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 90.0/len(sub_mask)
                freq = arange(40.+width/2,130.,width)
                if len(freq)>len(sub_mask):
                    freq = freq[0:-2]
                short_mask.append(sub_mask)
                short_mtime.append(time)        

short_data = array(short_data)
short_mask = array(short_mask)

if len(short_mask)==0:
    short_mask = zeros((len(short_data),len(short_data[0])))
    short_mtime = short_time

sortind = argsort(short_time)
sorttime = zeros(len(short_data))
sortshort = zeros((len(short_data),len(short_data[0])))
sortindm = argsort(short_mtime)
sortmask = zeros((len(short_data),len(short_data[0])))
for i in range(0,len(sortind)):
    sorttime[i] = short_time[sortind[i]]
    sortshort[i] = short_data[sortind[i]]
    sortmask[i] = short_mask[sortindm[i]]

mean_short, mean_mask = cf.time_mean(sortshort,sortmask)
mean_mask = ff.spike_flag(mean_short,mean_mask,freq,2.)
print freq[where(mean_mask==1.0)[0]]

short_array = ma.array(mean_short,mask=mean_mask)
short_comp = ma.compressed(short_array)
freq_array = ma.array(freq,mask=mean_mask)
freq_comp = ma.compressed(freq_array)

#Need to limit the frequencies for fit just to what we are considering.
fmin = where(freq_comp<=50.)[0][-1]
fmax = where(freq_comp<=100.)[0][-1]
print fmin, fmax
(Fa,Fb,Fc) = polyfit(freq_comp[fmin:fmax],short_comp[fmin:fmax],2)

savetxt(outdir+'June_'+date_ind+'_avg_short.txt',mean_short,delimiter=' ')
savetxt(outdir+'June_'+date_ind+'_fit_short.txt',polyval([Fa,Fb,Fc],freq),delimiter=' ')

#Plotting Checks
pylab.imshow(sortshort*10**9,aspect=90./len(sortind),extent=(40,130,len(sortind),0.0))
pylab.colorbar()
#pylab.title('Variation of Short over the day')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (short index)')
pylab.savefig(outdir+'June_'+date_ind+'_short_variation',dpi=300)
pylab.clf()

pylab.scatter(freq,mean_short*10**9,c='b',edgecolor='b',label='Unfiltered Mean Data')
pylab.plot(freq_comp,short_comp*10**9,label='Mean Data')
pylab.plot(freq,polyval([Fa,Fb,Fc],freq)*10**9,label='Linear Fit')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Power (nW)')
pylab.xlim(40,130)
pylab.legend()
pylab.grid()
#pylab.title('Mean Short and Fit')
pylab.savefig(outdir+'June_'+date_ind+'_short_mean',dpi=300)
pylab.clf()
