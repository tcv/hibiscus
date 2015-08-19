"""
Module to create gsm data for a given location and time. 
Set up to run in parallel using gsm_header.py and a batch script. 
Example batch script in /home/tabithav/scripts/mk_gsm.m on hippo.


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
import time

def main(gsm_times,antenna,fplot,site_params):
    start = time.time()
    antenna = int(antenna)

#Main directories and files for the input and output
    supdir = '../../supplemental_data/'
    if antenna==70:
        outdir = supdir+'gsm_data_70_Karoo_test/'
    elif antenna==100:
        outdir = supdir+'gsm_data_100_Karoo_test/'

#Parameters for the site.
    idate = site_params[0]
    lon = site_params[1]
    lat = site_params[2]
    elevation = site_params[3]
    gsm_freq = arange(50,111,1)

#Set up ra/dec arrays for GSM data
    gsm_var = numpy.load(supdir+'angles.npy')

#Set up az/alt arrays for antenna simulation data
    curr_az = numpy.load(supdir+'sim_az.npy')
    curr_alt = numpy.load(supdir+'sim_alt.npy')

#Load map data array
    gsm_data = numpy.load(supdir+'gsm_array.npy')

#Load beam simulation data array
    if antenna==70:
        sim_data = numpy.load(supdir+'beam_simulations50-90.npy')
    elif antenna==100:
        sim_data = numpy.load(supdir+'beam_simulations80-110.npy')
   
    print 'Initial loading time is: ',time.time()-start,' seconds'

    results = zeros(len(gsm_freq))
    find = 0

    tstart = time.time()

#Convert time to sidereal time. 
    origtime = float(gsm_times)/60.
    single_sid = ff.sidereal(origtime,idate,lon,lat)
    min = (single_sid%1.0)*60
    sec = (min%1.0)*60
    time_label = str(single_sid).split('.')[0]+'-'+str(min).split('.')[0]+'-'+str(sec).split('.')[0]
        
#Create sim ra/dec array for a given time. 
    sim_var = gf.radec_sim(curr_az,curr_alt,lat,lon,elevation,origtime,idate)

#Iterates over frequencies for a given time.        
    for f in range(0,len(gsm_freq)):
        fstart = time.time()
        freq = gsm_freq[f]

#Calculate expected Temperature for the given freq and time
        plot_label = outdir+'plots/gsm_sim_plots_'+time_label+'_hrs_'+str(int(freq))+'_MHz.png'
        fr = gf.ant_beam(gsm_data[f],gsm_var, sim_data[f], sim_var,plot_label,freq,float(fplot))
#            print 'Temperature value is: ', fr
#            print 'Frequency is: ', freq
        results[find] = fr
        find = find+1
#            print 'Running ',freq,' MHz took ',time.time()-fstart,' seconds.'

#Make a .npy array for the data
    if antenna==70:
        numpy.save(outdir+'gsm_data_Karoo_'+time_label+'_sid_time.npy',results)
    elif antenna==100: 
        numpy.save(outdir+'gsm_data_Karoo_'+time_label+'_sid_time.npy',results)
    find=0
#    print 'Running a single time sample takes: ',time.time()-tstart,' seconds'

    return

