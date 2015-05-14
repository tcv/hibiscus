"""
Module to create gsm data for a given location. 
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy
import pylab
#import scipy.interpolate as itp
#import numpy.ma as ma
#import scipy.optimize as opt
import os
#import skrf as rf
import sys
#import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('../../hibiscus'))
#import file_funcs as ff
#import gsm_funcs as gf


outdir = '../../supplemental_data/'
indir = '../../supplemental_data/gsm_data_70_Karoo/'

idate = '2015/4/1'
lon = '21.4109'
lat = '-30.7216'
elevation = 1080
gsm_freq = arange(50,110,1)

dirlist = os.listdir(indir)
full_data = zeros((len(dirlist),len(gsm_freq)))
full_times = zeros(len(dirlist))
tind= 0

for fname in dirlist:
    time_label = fname.split('_')[-3]
    time = float(time_label.split('-')[0])+float(time_label.split('-')[1])/60.+float(time_label.split('-')[2])/3600.
    data = numpy.load(indir+fname)
    full_data[tind] = data
    full_times[tind] = time
    tind = tind+1

sortind = argsort(full_times)
sorttime = zeros(len(full_times))
sortdata = zeros((len(full_times),len(gsm_freq)))
for i in range(0,len(full_times)):
    sorttime[i] = full_times[sortind[i]]
    sortdata[i] = full_data[sortind[i]]

pylab.imshow(sortdata,vmin=0,vmax = 5000,aspect = 060./24., extent=(gsm_freq[0],gsm_freq[-1],sorttime[-1],sorttime[0]))
pylab.colorbar()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.title('Expected Temperature (Kelvin)')
pylab.savefig(outdir+'GSM_Temp_70_Karoo.png')
pylab.clf()

pylab.plot(gsm_freq,sortdata[0])
pylab.plot(gsm_freq,sortdata[len(sortdata)/4])
pylab.plot(gsm_freq,sortdata[len(sortdata)/2])
pylab.plot(gsm_freq,sortdata[3*len(sortdata)/4])
pylab.xlim(50,110)
pylab.ylim(0,7000)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.savefig(outdir+'Sample_GSM_Temp_70_Karoo.png')
pylab.clf()

numpy.save(outdir+'gsm_data_full_70_Karoo.npy',sortdata)
numpy.save(outdir+'gsm_sid_time_full_70_Karoo.npy',sorttime)

