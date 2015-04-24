"""
Module for converting one hour of data into linear units,
truncating it to a set frequency range
and placing it into a single .npy array
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy
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
indir = '../../Karoo_data_Apr14_70/'
outdir = '../../Karoo_data_Apr14_70/'
directories = os.listdir(indir)

# Based upon the naming convention for the subdirectories in the raw data
for direct in directories:
#    direct = 'June'+date_ind+'_day_to_night'
    print direct
    dirlist = []
    date = direct.split('/')[-1]
    print date
    if len(direct.split('-'))>2:
        directory = indir+direct+'/'
        new_directory = directory
        dirlist = os.listdir(directory)
        hourdata = []
        noise = []
        short = []
        term = []
        load = []
        hourtime = []
        lt = []
        nt = []
        st = []
        tt = []
#Iterate for each file in the directory
        for fname in dirlist:
            if len(fname.split('_'))>3:
                filename = directory+fname
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
                mean_data = ma.mean(ma.compressed(ma.array(new_data,mask=new_mask)))

#Check for nan/inf (should be nulls)
                nandata = where(isnan(new_data))[0]
                for i in range(0,len(nandata)):
                    new_data[nandata[i]] = 0.0
                infdata = where(isinf(new_data))[0]
                for i in range(0,len(infdata)):
                    new_data[infdata[i]] = 0.0

                print 'Test Data value is: ',mean_data
                if fname.split('_')[-1]=='antenna.dat':
                    if mean_data> 1e-8 :
                        if mean_data<1e-7:
                            hourdata.append(new_data)
                            hourtime.append(time)
                elif fname.split('_')[-1]=='open.dat':
                    term.append(new_data)
                    tt.append(time)
                elif fname.split('_')[-1]=='short.dat':
                    short.append(new_data)
                    st.append(time)
                elif fname.split('_')[-1]=='50ohm.dat':
                    load.append(new_data)
                    lt.append(time)
                elif fname.split('_')[-1]=='noise.dat':
                    noise.append(new_data)
                    nt.append(time)
        numpy.save(new_directory+date+'_antenna.npy',hourdata)
        numpy.save(new_directory+date+'_ant_time.npy',hourtime)
        numpy.save(new_directory+date+'_term.npy',term)
        numpy.save(new_directory+date+'_term_time.npy',tt)
        numpy.save(new_directory+date+'_load.npy',load)
        numpy.save(new_directory+date+'_load_time.npy',lt)
        numpy.save(new_directory+date+'_short.npy',short)
        numpy.save(new_directory+date+'_short_time.npy',st)
        numpy.save(new_directory+date+'_noise.npy',noise)
        numpy.save(new_directory+date+'_noise_time.npy',nt)
        
