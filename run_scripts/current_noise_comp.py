"""
Module to calculate the 100 Ohm spectrum for a single day of data
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
Kdir='/lustre/tcv/calibration_data/'
directories = os.listdir(indir)

#Setting a single day for parallel computation
date_ind = sys.argv[1]

load_data = []
load_mask = []
load_time = []
load_mtime = []
term_data = []
term_mask = []
term_time = []
term_mtime = []
if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being rebinned is:',direct
    directory = indir+direct+'/'
    new_directory = outdir+direct+'/'
    dirlist = os.listdir(directory)
    short_data =loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')
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
            load_time.append(time) 
        elif fname.split('_')[-1]=='mask.dat':
            if fname.split('_')[-2]=='50ohm':
                filename = directory+fname
#load mask file 
                time,form,sub_mask,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 90.0/len(sub_mask)
                freq = arange(40.+width/2,130.,width)
                if len(freq)>len(sub_mask):
                    freq = freq[0:-2]
                load_mask.append(sub_mask)
                load_mtime.append(time)
            elif fname.split('_')[-2]=='open':
                filename = directory+fname
#load mask file
                time,form,sub_mask,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 90.0/len(sub_mask)
                freq = arange(40.+width/2,130.,width)
                if len(freq)>len(sub_mask):
                    freq = freq[0:-2]
                term_mask.append(sub_mask)
                term_mtime.append(time)

        elif fname.split('_')[-1]=='open.dat':
            filename = directory+fname
#load data file
            time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40.+width/2,130.,width)
            if len(freq)>len(sub_data):
                freq = freq[0:-2]
            term_data.append(sub_data)
            term_time.append(time)
        

load_data = array(load_data)
load_mask = array(load_mask)
term_data = array(term_data)
term_mask = array(term_mask)

sortind = argsort(load_time)
sorttime = zeros(len(load_data))
sortload = zeros((len(load_data),len(load_data[0])))
sortindm = argsort(load_mtime)
sortmask = zeros((len(load_data),len(load_data[0])))
for i in range(0,len(sortind)):
    sorttime[i] = load_time[sortind[i]]
    sortload[i] = load_data[sortind[i]]
    sortmask[i] = load_mask[sortindm[i]]

sortindt = argsort(term_time) 
sorttimet = zeros(len(term_data))
sortterm = zeros((len(term_data),len(term_data[0])))
sortindmt = argsort(term_mtime)
sortmaskt = zeros((len(term_data),len(term_data[0])))
for i in range(0,len(sortindt)):
    sorttimet[i] = term_time[sortindt[i]]
    sortterm[i] = term_data[sortindt[i]] 
    sortmaskt[i] = term_mask[sortindmt[i]]    

mean_load, mean_mask = cf.time_mean(sortload,sortmask)
mean_mask = ff.spike_flag(mean_load,mean_mask,freq,10.)
print "Bad 50 Ohm Frequencies are:", freq[where(mean_mask==1.0)[0]]

mean_term, mean_maskt = cf.time_mean(sortterm,sortmaskt)
mean_maskt = ff.spike_flag(mean_term,mean_maskt,freq,10.)
print "Bad 100 Ohm Frequencies are:", freq[where(mean_maskt==1.0)[0]] 

cal_mean_load = mean_load*300./(mean_load-short_data)
cal_mean_term = mean_term*300./(mean_load-short_data) 
cal_mean_short = short_data*300./(mean_load-short_data)

load_array = ma.array(cal_mean_load,mask=mean_mask)
load_comp = ma.compressed(load_array)
freq_array = ma.array(freq,mask=mean_mask)
freq_comp = ma.compressed(freq_array)
term_array = ma.array(cal_mean_term,mask=mean_maskt)
term_comp = ma.compressed(term_array)
freq_arrayt = ma.array(freq,mask=mean_maskt)
freq_compt = ma.compressed(freq_arrayt)
short_array = ma.array(cal_mean_short,mask=mean_mask)
short_comp = ma.compressed(short_array)

#pylab.scatter(freq,cal_mean_load,c='b',edgecolor='b',s=5,label='Unfiltered Mean Data')
pylab.plot(freq_comp,load_comp,'b.-',label='50 Ohm')
#pylab.scatter(freq,cal_mean_term,c='g',edgecolor='g',s=5,label='Unfiltered Mean Data')
pylab.plot(freq_compt,term_comp,'g.-',label='100 Ohm')
#pylab.plot(freq,polyval([Fa,Fb,Fc],freq)*10**9,'g',label='Linear Fit')
pylab.plot(freq_comp,short_comp,'r.-',label='Short')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.xlim(60,100)
pylab.legend()
pylab.grid()
#pylab.title('Mean 100 Ohm and Fit')
pylab.savefig(outdir+'June_'+date_ind+'_compare_mean',dpi=300)
pylab.clf()
