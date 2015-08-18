"""
Module to make plots for each hour dataset.
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


gbt_idate = '2013/5/1'
gbt_lat = '38.433'
gbt_lon = '-79.8397'

guad_idate = '2013/6/1'
guad_lat = '28.8833'
guad_lon = '-118.3'

karoo_idate = '2015/4/1'
karoo_lat = '-30.7216'
karoo_lon = '21.4109'

#indir = '../../Karoo_data_npy/'
#plotdir = '../../Karoo_data_npy_plts/'
indir = sys.argv[1]
plotdir = sys.argv[2]
#plotdir=indir
directories = os.listdir(indir)
total_average = []
total_dates = []

for direct in directories:
    if direct.split('_')[-1]=='100':
        dirlist = os.listdir(indir+direct+'/')
        sub_dir = indir+direct+'/'
        plot_sub = plotdir+direct+'/'
        print plot_sub
        ff.mkdir_p(str(plot_sub))
        day = direct.split('_')[0]
        total_dates.append(day)
        average_data = []
        for fname in dirlist:
            single = []
            mean = []
            if fname.split('_')[-1]=='antenna.npy':
                date = fname.split('_')[0]
                print date            
                single = numpy.load(sub_dir+fname)
                if len(single)<2000:
                    mean = ma.mean(single,axis=0)*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    pylab.imshow(single*1e9,vmin=0,vmax=200,aspect = 90./len(single),extent=(40,130,len(single),0))
                    pylab.colorbar()
                    pylab.xlabel('Frequency (MHz)')
                    pylab.ylabel('Time (sample)')
                    pylab.savefig(plot_sub+date+'_antenna_waterfall.png',dpi=300)
                    pylab.clf()
                    num_of_datasets = str(len(single))
                else:
                    mean = single*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    num_of_datasets = str(1)
                pylab.plot(freqs,mean)
                pylab.ylim(0,200)
                pylab.xlabel('Frequency (MHz)')
                pylab.ylabel('Power (nW)')
                pylab.grid()
                pylab.title('Mean of '+num_of_datasets+' Datasets')
                pylab.savefig(plot_sub+date+'_antenna_mean.png',dpi=300)
                pylab.clf()
                average_data.append(mean)
            elif fname.split('_')[-1]=='short.npy':
                date = fname.split('_')[0]
                print date            
                single = numpy.load(sub_dir+fname)
                if len(single)<2000:
                    mean = ma.mean(single,axis=0)*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    pylab.imshow(single*1e9,vmin=0,vmax=100,aspect = 90./len(single),extent=(40,130,len(single),0))
                    pylab.colorbar()
                    pylab.xlabel('Frequency (MHz)')
                    pylab.ylabel('Time (sample)')
                    pylab.savefig(plot_sub+date+'_short_waterfall.png',dpi=300)
                    pylab.clf()
                    num_of_datasets = str(len(single))
                else:
                    mean = single*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    num_of_datasets = str(1)
                pylab.plot(freqs,mean)
                pylab.ylim(0,100)
                pylab.xlabel('Frequency (MHz)')
                pylab.ylabel('Power (nW)')
                pylab.grid()
                pylab.title('Mean of '+num_of_datasets+' Datasets')
                pylab.savefig(plot_sub+date+'_short_mean.png',dpi=300)
                pylab.clf()
            elif fname.split('_')[-1]=='load.npy':
                date = fname.split('_')[0]
                print date            
                single = numpy.load(sub_dir+fname)
                if len(single)<2000:
                    mean = ma.mean(single,axis=0)*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    pylab.imshow(single*1e9,vmin=0,vmax=100,aspect = 90./len(single),extent=(40,130,len(single),0))
                    pylab.colorbar()
                    pylab.xlabel('Frequency (MHz)')
                    pylab.ylabel('Time (sample)')
                    pylab.savefig(plot_sub+date+'_load_waterfall.png',dpi=300)
                    pylab.clf()
                    num_of_datasets = str(len(single))
                else:
                    mean = single*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    num_of_datasets = str(1)
                pylab.plot(freqs,mean)
                pylab.ylim(0,100)
                pylab.xlabel('Frequency (MHz)')
                pylab.ylabel('Power (nW)')
                pylab.grid()
                pylab.title('Mean of '+num_of_datasets+' Datasets')
                pylab.savefig(plot_sub+date+'_load_mean.png',dpi=300)
                pylab.clf()
            elif fname.split('_')[-1]=='term.npy':
                date = fname.split('_')[0]
                print date            
                single = numpy.load(sub_dir+fname)
                if len(single)<2000:
                    mean = ma.mean(single,axis=0)*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    pylab.imshow(single*1e9,vmin=0,vmax=100,aspect = 90./len(single),extent=(40,130,len(single),0))
                    pylab.colorbar()
                    pylab.xlabel('Frequency (MHz)')
                    pylab.ylabel('Time (sample)')
                    pylab.savefig(plot_sub+date+'_term_waterfall.png',dpi=300)
                    pylab.clf()
                    num_of_datasets = str(len(single))
                else:
                    mean = single*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    num_of_datasets = str(1)
                pylab.plot(freqs,mean)
                pylab.ylim(0,100)
                pylab.xlabel('Frequency (MHz)')
                pylab.ylabel('Power (nW)')
                pylab.grid()
                pylab.title('Mean of '+num_of_datasets+' Datasets')
                pylab.savefig(plot_sub+date+'_term_mean.png',dpi=300)
                pylab.clf()
            elif fname.split('_')[-1]=='noise.npy':
                date = fname.split('_')[0]
                print date            
                single = numpy.load(sub_dir+fname)
                if len(single)<2000:
                    mean = ma.mean(single,axis=0)*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    pylab.imshow(single*1e9,vmin=0,vmax=100,aspect = 90./len(single),extent=(40,130,len(single),0))
                    pylab.colorbar()
                    pylab.xlabel('Frequency (MHz)')
                    pylab.ylabel('Time (sample)')
                    pylab.savefig(plot_sub+date+'_noise_waterfall.png',dpi=300)
                    pylab.clf()
                    num_of_datasets = str(len(single))
                else:
                    mean = single*1e9
                    freqs = arange(40.,130.,90./len(mean))
                    num_of_datasets = str(1)
                pylab.plot(freqs,mean)
                pylab.ylim(0,100)
                pylab.xlabel('Frequency (MHz)')
                pylab.ylabel('Power (nW)')
                pylab.grid()
                pylab.title('Mean of '+num_of_datasets+' Datasets')
                pylab.savefig(plot_sub+date+'_noise_mean.png',dpi=300)
                pylab.clf()

        total_average.append(average_data)

#numpy.save(indir+'total_average.npy',total_average)
#numpy.save(indir+'total_dates.npy',total_dates)

for i in range(0,len(total_average)):
    day_data = total_average[i]
    for j in range(0,len(day_data)):
        pylab.plot(freqs,day_data[j])
    pylab.ylim(0,200)
    pylab.ylabel('Power (nW)')
    pylab.xlabel('Frequency (MHz)')
    pylab.grid()
    pylab.savefig(plotdir+total_dates[i]+'_antenna_averages.png',dpi=300)
    pylab.clf()
    pylab.imshow(day_data,vmin=0,vmax=200,aspect = 90./len(day_data),extent=(40,130,len(day_data),0))
    pylab.colorbar()
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Time (hour blocks)')
    pylab.savefig(plotdir+total_dates[i]+'_antenna_avg_waterfall.png',dpi=300)
    pylab.clf()
    
