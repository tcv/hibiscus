"""
Module for truncating the data in frequency and converting to linear units.
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import sys
sys.path.append(os.path.abspath('../../hibiscus'))
import file_funcs as ff

# Main directories for the input and output
#indir = '/lustre/anat/Guadalupe_data/'
indir = '../../Karoo_data_Apr14_70/'
outdir = '../../Karoo_data_Apr14_70/truncated/'
#outdir = '/lustre/tcv/truncated_data/'
directories = os.listdir(indir)

#Setting a single day for parallel computation
#date_ind = sys.argv[1]

#if int(date_ind)<15:
# Based upon the naming convention for the subdirectories in the raw data
for direct in directories:
#    direct = 'June'+date_ind+'_day_to_night'
    print direct
    dirlist = []
    if len(direct.split('-'))>2:
        directory = indir+direct+'/'
        new_directory = outdir
        dirlist = os.listdir(directory)
#Iterate for each file in the directory
    for fname in dirlist:
        if len(fname.split('_'))>=3:
            filename = directory+fname
#            print filename
            if os.path.isfile(indir+'truncated/'+fname)==False:
#load data file
                print filename
                time,form,sub_data,mask,freq,volt,temp = ff.loadsingle(filename)
                width = 250.0/len(sub_data)
                freqs = arange(0,250.0,width)
                mask = zeros(len(sub_data))
#truncate to limited frequencies
                new_data,new_mask,new_freq = ff.truncate(sub_data,mask,freqs,40.,130.)
                new_data = array(new_data)
#Convert to linear units
                new_data = 10**(new_data/10.)

#Check for nan/inf (should be nulls)
                nandata = where(isnan(new_data))[0]
                for i in range(0,len(nandata)):
                    new_data[nandata[i]] = 0.0
                infdata = where(isinf(new_data))[0]
                for i in range(0,len(infdata)):
                    new_data[infdata[i]] = 0.0

#Output linear data within voltage limits 
#            if volt>10.0:
#                if volt<12.0:
                ff.writesingle(filename,new_directory,new_data,'')

#For if there is no voltage information.
#            elif volt==0.0:
#                ff.writesingle(filename,new_directory,new_data,'')
            
