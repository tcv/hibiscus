"""
Module to rebin a single day of data
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
sys.path.append(os.path.abspath('../../hibiscus'))
import file_funcs as ff

#Main directories for the input and output
#indir = '/lustre/tcv/truncated_data/'
#indir = '../../Reverb_chamber_data/'
indir = '../../Karoo_data_day1/truncated/'
#outdir = '/lustre/tcv/freq_rebinned_data/'
outdir = '../../Karoo_data_day1/rebinned/'
#Kdir = '/lustre/tcv/rfi_check_data/'
#directories = os.listdir(indir)

#Setting a single day for parallel computation
#date_ind = sys.argv[1]

#Setting Rebinning scales
timescale = 32
freqscale = 32

total_mask = 0.
mid_mask = 0.
sub_mid_mask = 0.
sub_tot_mask = 0.
total_sum = 0.
total_sub_sum = 0.

#if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
#    direct = 'June'+date_ind+'_day_to_night'
#    print 'Directory being rebinned is:',direct
#    directory = indir+direct+'/'
#    new_directory = outdir+direct+'/'
new_directory = outdir
#if freqscale==32:
#repping calibration data for use
#    short_data = loadtxt(Kdir+'June_'+date_ind+'_fit_short.txt')
#    short_full = loadtxt(Kdir+'June_'+date_ind+'_avg_short.txt')

#    dirlist = os.listdir(directory)
dirlist = os.listdir(indir)
    #print dirlist
#Iterate for each file in the directory
for fname in dirlist:
    if len(fname.split('-'))>=5:
        if fname.split('-')[0]=='2015':
#            if fname.split('-')[1]=='.dat':
#            if fname.split('-')[-1]!='cal.dat':
#                filename = directory+fname
            filename = indir + fname
            if os.path.isfile(outdir+fname)==False:
                print fname
#load data file
                time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 90.0/len(sub_data)
                freq = arange(40,130.0,width)
#                freq = arange(40,130.0,width)
#                fmin = where(freq<=50.)[0][-1]
#                fmax = where(freq<=100.)[0][-1]

#basic freq flagging
                mask = ff.flagging(sub_data,freq,3.,freqscale)
                spike_mask = ff.spike_flag(sub_data,mask,freq,100.)
#                short_mask = ff.cal_flag(short_full,short_data,spike_mask,freq,1.e-10)
#                thresh_mask = ff.threshold_flag(sub_data,short_mask,freq,5.)
#                mid_mask = mid_mask+sum(mask)
#                sub_mid_mask = sub_mid_mask+sum(mask[fmin:fmax])
#                for i in range(0,len(sub_data)):
#                    if spike_mask[i]==1.0:
#                        mask[i] = 1.0
#                total_mask = total_mask+sum(spike_mask)
#                mask = thresh_mask
                mask = spike_mask
                total_mask = total_mask+sum(mask)
#                sub_tot_mask = sub_tot_mask+sum(mask[fmin:fmax])
                total_sum = total_sum +len(mask)
#                total_sub_sum = total_sub_sum+len(mask[fmin:fmax])
                
#freq rebinning
                new_data,new_mask,new_freq = ff.rebin(sub_data,mask,freq,freqscale) 
                new_data = array(new_data)
                new_mask = array(new_mask)

#Check for nan/inf (should be nulls)
                nandata = where(isnan(new_data))[0]
                for i in range(0,len(nandata)):
                    new_data[nandata[i]] = 0.0
                infdata = where(isinf(new_data))[0]
                for i in range(0,len(infdata)):
                    new_data[infdata[i]] = 0.0
 
                ff.writesingle(filename,new_directory,new_data,'')
                ff.writesingle(filename,new_directory,new_mask,'_mask')

#print 'Percent of Data Flagged:',100.*total_mask/total_sum
print 'Percent of Data Flagged with basic flagger: ',100.*mid_mask/total_sum
print 'Percent of Data Flagged with additional spike flaggers: ',100*total_mask/total_sum
#print 'Percent of Data between 50-100 MHz with basic flagger: ',100.*sub_mid_mask/total_sub_sum
#print 'Percent of Data between 50-100 MHz with add flagger: ',100.*sub_tot_mask/total_sub_sum

