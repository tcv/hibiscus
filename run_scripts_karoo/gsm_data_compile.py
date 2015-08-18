"""
Module to combine gsm data for all sidereal times into a single file for a given location and antenna.

Also makes plots of the gsm model data. 

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


#outdir = '../../supplemental_data/'
indir = '../../supplemental_data/gsm_data_100_Karoo_test/'
outdir = indir
idate = '2015/4/1'
lon = '21.4109'
lat = '-30.7216'
elevation = 1080
gsm_freq = arange(50,111,1)

dirlist = os.listdir(indir)
full_data = zeros((len(dirlist),len(gsm_freq)))
full_times = zeros(len(dirlist))
tind= 0

for fname in dirlist:
    if fname.split('.')[-1]=='npy':
        if fname.split('_')[2]=='Karoo':
            time_label = fname.split('_')[-3]
            time = float(time_label.split('-')[0])+float(time_label.split('-')[1])/60.+float(time_label.split('-')[2])/3600.
            data = numpy.load(indir+fname)
            full_data[tind] = data
            full_times[tind] = time
            tind = tind+1
full_data = full_data[0:tind]
full_times = full_times[0:tind]
#full_data = numpy.load(indir+'gsm_data_full_100_Karoo.npy')
#full_times = numpy.load(indir+'gsm_sid_time_full_100_Karoo.npy')
print min(full_times),max(full_times)

sortind = argsort(full_times)
sorttime = zeros(len(full_times))
sortdata = zeros((len(full_times),len(gsm_freq)))
for i in range(0,len(full_times)):
    sorttime[i] = full_times[sortind[i]]
    sortdata[i] = full_data[sortind[i]]

#sortdata = numpy.load(indir+'gsm_data_full_70_Karoo.npy')
#sorttime = numpy.load(indir+'gsm_sid_time_full_70_Karoo.npy')


pylab.imshow(sortdata[:,30:-1],vmin=0,vmax = 10000,aspect = 030./24., extent=(gsm_freq[30],gsm_freq[-1],sorttime[-1],sorttime[0]))
pylab.colorbar()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.title('Expected Temperature (Kelvin)')
pylab.savefig(outdir+'GSM_Temp_100_Karoo.png')
pylab.clf()

pylab.plot(gsm_freq,sortdata[0])
pylab.plot(gsm_freq,sortdata[len(sortdata)/4])
pylab.plot(gsm_freq,sortdata[len(sortdata)/2])
pylab.plot(gsm_freq,sortdata[3*len(sortdata)/4])
pylab.xlim(80,110)
pylab.ylim(0,4000)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.savefig(outdir+'Sample_GSM_Temp_100_Karoo.png')
pylab.clf()


pylab.plot(sorttime,sortdata[:,44],label='95 MHz')
pylab.plot(sorttime,sortdata[:,54],label='105 MHz')
pylab.plot(sorttime,sortdata[:,34],label='85 MHz')
pylab.xlim(0,24.)
pylab.xlabel('Sidereal Time (Hours)')
pylab.ylim(0,5e3)
pylab.ylabel('Temperature (Kelvin)')
pylab.grid()
pylab.legend()
pylab.savefig(outdir+'Sample_GSM_Time_stream_100_Karoo.png')
pylab.clf()

numpy.save (outdir + 'gsm_data_full_100_Karoo.npy', sortdata)
numpy.save(outdir+'gsm_sid_time_full_100_Karoo.npy',sorttime)

