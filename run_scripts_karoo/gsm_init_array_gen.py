"""
Code for making the initial loading files for gsm_generate_opt.py
Converts .dat files data into .npy arrays
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

def mk_gsm_files(indir,outdir):
    start = time.time()
#Main directories and files for the input and output
    supdir = indir

    radecfile = supdir+'angles.dat'
    gsmdir = supdir+'galaxy_maps_radec/'

    gsm_freq = arange(50,111,1)

#Set up ra/dec arrays for GSM data
    ras,decs = gf.get_gsm_radec(radecfile)
    gsm_var= vstack((ras,decs)).T
    numpy.save(outdir+'angles.npy',gsm_var)

#Set up arrays for gsm data:
    gsm_data = zeros((len(gsm_freq),len(ras)))

    for i in range(0,len(gsm_freq)):
        freq = gsm_freq[i]
        gsm_ind = arange(0,len(ras))
        gsmdata = gf.gsm_temps(gsmdir,freq,gsm_ind)
        gsm_data[i] = gsmdata

    numpy.save(outdir+'gsm_array.npy',gsm_data)

    print 'GSM file creation time is: ',time.time()-start,' seconds'

    return

def mk_sim_files(antenna,indir,outdir):
    start = time.time()
    antenna = int(antenna)
#Main directories and files for the input and output
    supdir = indir
    beamfile = supdir+'beam_simulations50-90.dat'

#Parameters for the site.
    idate = '2015/4/1'
    lon = '21.4109'
    lat = '-30.7216'
    elevation = 1080
    gsm_freq = arange(50,111,1)
 
#Set up arrays for sim data:
    sim_data = -50.*ones((len(gsm_freq),181*181*2))
 
    for i in range(0,len(gsm_freq)):
        freq = gsm_freq[i]
#Get HIbiscus beam sim data (dB), theta and phi for a given frequency (theta is alt, phi is az).
        curr_az,curr_alt,gaindb = gf.antenna_beam_full(beamfile,freq)
        if antenna==70:
            sim_data[i] = gaindb
        elif antenna==100:
            if i==0:
                lin_fit = (gaindb-sim_data[0])/(gsm_freq[30]-gsm_freq[0])
                sim_data[i] = lin_fit*gsm_freq[i]+(sim_data[0]-lin_fit*gsm_freq[0])
                sim_data[i+30] = gaindb
            elif freq<=80.:
                lin_fit = (sim_data[30]-sim_data[0])/(gsm_freq[30]-gsm_freq[0])
                sim_data[i] = lin_fit*gsm_freq[i]+sim_data[0]-lin_fit*gsm_freq[0]
                sim_data[i+30] = gaindb
 
    numpy.save(outdir+'sim_az.npy',curr_az)
    numpy.save(outdir+'sim_alt.npy',curr_alt)
    if antenna==70:
        numpy.save(outdir+'beam_simulatons50-90.npy',sim_data)
    elif antenna==100:
        numpy.save(outdir+'beam_simulations80-110.npy',sim_data)
 
    print 'Beam Data Array Creation time is: ',time.time()-start,' seconds'
 
    return


antenna = '100'
supdir = '../../supplemental_data/'
outdir = supdir
start = time.time()
mk_gsm_files(supdir,outdir)
mk_sim_files(antenna,supdir,outdir)
print 'Code took ',time.time()-start,' seconds to run.'


