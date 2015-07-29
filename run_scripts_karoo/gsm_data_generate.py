"""
Module to create gsm data for a given location and time. 
Set up to run in parallel using python_for_bash.py and a batch script. 
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

def main(gsm_times):

#Main directories for the input and output
    outdir = '../../supplemental_data/gsm_data_70_Karoo/'
    beamfile = '../../supplemental_data/beam_simulations50-90.dat'
    radecfile = '../../supplemental_data/angles.dat'
    gsmdir = '../../supplemental_data/galaxy_maps_radec/'

#Parameters for the site.
    idate = '2015/4/1'
    lon = '21.4109'
    lat = '-30.7216'
    elevation = 1080
    gsm_freq = arange(50,111,1)

# Will need to add an offset due to phi=0 at magnetic south
# Rather than geographic south (like currently assumed). 
   
#Convert time to sidereal time. 
    time = float(gsm_times)/60.
    single_sid = ff.sidereal(time,idate,lon,lat)
    print 'Sidereal Time is: ' ,single_sid
    min = (single_sid%1.0)*60 
    sec = (min%1.0)*60
    time_label = str(single_sid).split('.')[0]+'-'+str(min).split('.')[0]+'-'+str(sec).split('.')[0]
    print time_label

#Set up ra/dec arrays for GSM data
    gsm_data = zeros(len(gsm_freq))
    ras,decs = gf.get_gsm_radec(radecfile)
    find = 0
    tind = 0

    for freq in gsm_freq:

#Get HIbiscus beam sim data (dB), theta and phi for a given frequency (theta is alt, phi is az).
        gaindb, sim_var = gf.sim_comp(beamfile,freq)

#Load GSM file for a given frequency
        gsm_ind = arange(0,len(ras))
        gsmdata = gf.gsm_temps(gsmdir,freq,gsm_ind)

#Get GSM, az, alt arrays for a given freq, sidereal time
        gsm_array, gsm_var = gf.gsm_comp(gsmdata,ras,decs,lat,lon,elevation,time,idate)
        print shape(gsm_array), shape(gsm_var)
        print shape(gaindb),shape(sim_var)

#Calculate expected Temperature for the given freq, time
        plot_label = outdir+'plots/gsm_sim_plots_'+time_label+'_hrs_'+str(int(freq))+'_MHz.png'
        fr = gf.ant_beam(gsm_array,gsm_var, gaindb, sim_var,plot_label,freq)
        print 'Temperature value is: ', fr
        print 'Frequency is: ', freq
        gsm_data[find] = fr
        find = find+1

#Make a .npy array for a given sidereal time 
    print shape(gsm_data)
    numpy.save(outdir+'gsm_data_Karoo_'+time_label+'_sid_time.npy',gsm_data)

    return

