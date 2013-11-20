"""
Module to rebin a single day of data
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
#from pylab import *
#import matplotlib.pyplot as plt
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import sys
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

maindir = '/lustre/tcv/truncated_data/'
outdir = '/lustre/tcv/freq_rebinned_data/'

timescale = 32
binscale = 32

directories = os.listdir(maindir)

date_ind = sys.argv[1]
total_mask = 0.
mid_mask = 0.
total_sum = 0.
if int(date_ind)<15:
    direct = 'June'+date_ind+'_day_to_night'
    print 'Directory being rebinned is:',direct
    directory = maindir+direct+'/'
    new_directory = outdir+direct+'/'
    dirlist = os.listdir(directory)
    for fname in dirlist:
        if fname.split('_')[-2] =='antenna':
            filename = directory+fname
            time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(filename)
            width = 90.0/len(sub_data)
            freq = arange(40,130.0,width)
            mask = fc.flagging(sub_data,freq,3.,binscale)
            spike_mask = fc.spike_flag(sub_data,100.)
            mid_mask = mid_mask+sum(mask)
            for i in range(0,len(sub_data)):
                if spike_mask[i]==1.0:
                    mask[i] = 1.0
            total_mask = total_mask+sum(mask)
            total_sum = total_sum +len(mask)
            new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freq,binscale)       
            new_data = array(new_data)
            nandata = where(isnan(new_data))[0]
            for i in range(0,len(nandata)):
                new_data[nandata[i]] = 0.0
            infdata = where(isinf(new_data))[0]
            for i in range(0,len(infdata)):
                new_data[infdata[i]] = 0.0
            fc.writesingle(filename,new_directory,new_data)

print 'Percent of Data Flagged:',100.*total_mask/total_sum
print 'Percent of Data Flagged without Spike Flagger:',100.*mid_mask/total_sum
