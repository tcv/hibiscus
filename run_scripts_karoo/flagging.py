"""
Module to generated mask for a given data set. 
Takes 4 inputs (directory for real data, earliest time file day and hour in that directory, and hour between 0 and 23.)
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
import os
import skrf as rf
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('../../hibiscus'))
import file_funcs as ff
import gsm_funcs as gf
import cal_funcs as cf
import eff_funcs as ef
import ephem as eph
import errno
import numpy.polynomial.polynomial as poly

indir = sys.argv[1]
outdir = indir
#outdir= '../../supplemental_data/'
directories = os.listdir(indir)

#Load files for each hour full day of data. 
#Set to max length, then truncate. Max length is 1 dataset every 3 seconds for 24 hours. 
#Freq length set by previous experience, and is 11796 (truncated, but not rebinned).
day_idate = int(sys.argv[2])*24+int(sys.argv[3])
fscale = 32
tscale = 32
data = zeros((1200,11796))
mask = zeros((1200,11796))
times = zeros((1200))
dint = 0
tint = 0

hour = int(sys.argv[4])
for direct in directories:
    curr_day =int(direct.split('-')[2])
    end = direct.split('-')[-1]
    cur_hour = int(end.split('_')[0])
    cur_date = curr_day*24+cur_hour
    if (cur_date == day_idate+hour):
        if direct.split('_')[-1]=='antenna.npy':
            single = numpy.load(indir+direct) 
            if len(single)<2000:
                freqs = arange(40.,130.,90./len(single[0]))
                for time in range(0,len(single)):
                    data[dint] = single[time]
                    first_mask = ff.flagging(single[time],freqs,3.,fscale)
                    spike_mask = ff.spike_flag(single[time],first_mask,freqs, 100.)
                    mask[dint] = spike_mask
                    dint+=1
             
                file =  str(direct.split('/')[-1])
                date = str(file.split('_')[0])
                print shape(single)
        elif direct.split('_')[-2] =='ant':
            single = numpy.load(indir+direct)
            for time in range(0,len(single)):
                times[tint] = single[time]
                tint+=1

data = data[0:dint]
mask = mask[0:dint]
times = times[0:tint]


for f in range(0,len(freqs)):
    new_mask = ff.timeflag(data[:,f],mask[:,f],times,3.,tscale)
    new_mask = ff.threshold_flag(data[:,f],new_mask[:,f],freqs[f],50.)
    mask[:,f] = new_mask

numpy.save(outdir+'/'+date+'_mask.npy',mask)
print shape(mask)
print 'Percent of Data Flagged from Frequency and Time Masking: ',100.*sum(mask)/(len(mask)*len(mask[0]))

