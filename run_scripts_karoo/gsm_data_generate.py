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

def main(gsm_times):
#Main directories for the input and output
    outdir = '../../supplemental_data/'
    beamfile = '../../supplemental_data/beam_simulations50-90.dat'
    radecfile = '../../supplemental_data/angles.dat'
    gsmdir = '../../supplemental_data/galaxy_maps_radec/'


    idate = '2015/4/1'
    lon = '21.4109'
    lat = '-30.7216'
    elevation = 1080
    gsm_freq = arange(50,111,1)
    
#gsm_freq = sys.argv[1]
#gsm_times = arange(0,24.,0.05)
#    gsm_times = sys.argv[1]
#print shape(gsm_times), shape(gsm_freq)

    time = float(gsm_times)/60.
    single_sid = ff.sidereal(time,idate,lon,lat)
    print 'Sidereal Time is: ' ,single_sid
    min = (single_sid%1.0)*60 
#print min 
    sec = (min%1.0)*60

    time_label = str(single_sid).split('.')[0]+'-'+str(min).split('.')[0]+'-'+str(sec).split('.')[0]
    print time_label

#gsm_data = zeros((len(gsm_freq),len(gsm_times)))
#gsm_data = zeros(len(gsm_times))
    gsm_data = zeros(len(gsm_freq))
#sid_times = zeros(len(gsm_times))
    ras,decs = gf.get_gsm_radec(radecfile)
#    nras,ndecs,gsm_ind = gf.trunc_gsm(ras,decs,9,18)
#    print shape(nras),shape(ndecs),shape(gsm_ind)
#time = float(gsm_times)
    find = 0
    tind = 0
#gsm_array, gsm_var = gf.gsm_comp(gsmdata,ras,decs,lat,lon,elevation,time,idate)
    for freq in gsm_freq:
        gaindb, sim_var = gf.sim_comp(beamfile,freq)
        gsm_ind = arange(0,len(ras))
        freq_gsm, gsmdata = gf.gsm_temps(gsmdir,freq,gsm_ind)
#for time in gsm_times:

        gsm_array, gsm_var = gf.gsm_comp(gsmdata,ras,decs,lat,lon,elevation,time,idate)


        print shape(gsm_array), shape(gsm_var)
        print shape(gaindb),shape(sim_var)
#        adj = 0.001
#        grid_alt, grid_az = numpy.mgrid[0+adj:pi/2.-adj:90j,0+adj:2.*pi-adj:180j]
#        print shape(grid_alt),shape(grid_az) 
#        print where(isnan(gsm_array))
#        print where(isnan(gsm_var))
#        print where(isinf(gsm_array))
#        print where(isinf(gsm_var))
#        print gsm_array, gsm_var
        print where(isnan(sim_var)),where(isnan(gaindb))
        print where(isinf(sim_var)),where(isnan(gaindb))

        fr = gf.ant_beam(gsm_array,gsm_var, gaindb, sim_var)
        print 'Temperature value is: ', fr
        print 'Frequency is: ', freq
#    gsm_data[find,tind] = fr
#    gsm_data[tind] = fr
        gsm_data[find] = fr
        find = find+1
#single_sid = ff.sidereal(time,idate,lon,lat)
#    sid_times[tind] = single_sid
#print 'Sidereal Time is: ' ,single_sid
#    tind = tind+1
#    find = find+1
#    tind = 0

    print shape(gsm_data)
#min = (single_sid%1.0)*60
#print min
#sec = (single_sid%1.0)*60

#time_label = str(single_sid).split('.')[0]+'-'+str(min).split('.')[0]+'-'+str(sec).split('.')[0]

#print time_label
    numpy.save(outdir+'gsm_data_Karoo_'+time_label+'_sid_time.npy',gsm_data)

#numpy.save(outdir+'gsm_data_Karoo_'+str(int(gms_freq))+ '_MHz.npy',gsm_data)
#numpy.save(outdir+'gsm_times_Karoo_'+str(int(gsm_freq))+'_MHz.npy',sid_times)

    return

