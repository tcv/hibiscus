"""
Module to create gsm data for a given location. 
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

#Main directories for the input and output
outdir = '../../Karoo_data_Apr12_70/'
beamfile = '../../beam_simulations50-90.dat'
radecfile = '../../angles.dat'
gsmdir = '../../galaxy_maps_radec/'



idate = '2015/4/1'
lon = '21.4109'
lat = '-30.7216'
elevation = 1080
gsm_freq = arange(50,110,1)
gsm_times = arange(0,24.,0.05)
print shape(gsm_times), shape(gsm_freq)


gsm_data = []
sid_times = []
for freq in gsm_freq:
    gaindb, sim_var = gf.sim_comp(beamfile,freq)
    freq_gsm, gsmdata = gf.gsm_temps(gsmdir,freq)
    ras,decs = gf.get_gsm_radec(radecfile)
    single_freq = []
    for time in gsm_times:
        gsm_array, gsm_var = gf.gsm_comp(gsmdata,ras,decs,lat,lon,elevation,time,idate)
        fr = gf.ant_beam(gsm_array,gsm_var, gaindb, sim_var)
        print fr
        single_freq.append(fr)
        if freq == gsm_freqs[0]:
            single_sid = ff.sidereal(time,idate,lon,lat)
            sid_times.append(single_sid)
    gsm_data.append(single_freq)

print shape(gsm_data)

numpy.save(outdir+'gsm_data_Karoo.npy',gsm_data)
numpy.save(outdir+'gsm_times_Karoo.npy',sid_times)



