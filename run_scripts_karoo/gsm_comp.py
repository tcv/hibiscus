"""
Module to do a rough comparison of gsm and real data.
Takes 3 inputs (directory for real data, earliest time file day and hour in that directory)
Creates a plot to show the comparison at 70 MHz.
Includes a best guess for the calibration factor at 70 MHz to make the plot work. 
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

indir = sys.argv[1]
outdir= '../../supplemental_data/'
supdir = '../../supplemental_data/'


gbt_idate = '2013/5/1'
gbt_lat = '38.433'
gbt_lon = '-79.8397'

guad_idate = '2013/6/1'
guad_lat = '28.8833'
guad_lon = '-118.3'

karoo_idate = '2015/4/1'
karoo_lat = '-30.7216'
karoo_lon = '21.4109'

gsm_data = numpy.load(supdir+'gsm_data_full_70_Karoo.npy')
gsm_times = numpy.load(supdir+'gsm_sid_time_full_70_Karoo.npy')
gsm_freqs = arange(50,111,1)
directories = os.listdir(indir)
data = []
times = []
day_idate = int(sys.argv[2])*24+int(sys.argv[3])

for hour in range(0,24):
    for direct in directories:
         curr_day =int(direct.split('-')[2])
         end = direct.split('-')[-1]
         cur_hour = int(end.split('_')[0])
         cur_date = curr_day*24+cur_hour
         if (cur_date == day_idate+hour):
             if direct.split('_')[-1]=='antenna.npy':
                 single = numpy.load(indir+direct) 
                 if len(single)<2000:
                     freqs = arange(40.,130.,90./len(single[0]))
                     for time in range(0,len(single)):
                         data.append(single[time])
             elif direct.split('_')[-2] =='ant':
                 single = numpy.load(indir+direct)
                 for time in range(0,len(single)):
                     times.append(single[time]-7.3)
print shape(data)
print shape(times)

sid_times = zeros(len(times))
for i in range(0,len(times)):
    single_sid = ff.sidereal(times[i],karoo_idate,karoo_lon,karoo_lat)
    sid_times[i] = single_sid

f70_gsm = where(gsm_freqs<=70.)[0][-1]
f70_data = where(freqs<=70.)[0][-1]
print f70_gsm, f70_data

data_70 = zeros(len(times))
for i in range(0,len(times)):
    data_70[i] = data[i][f70_data]

pylab.plot(gsm_times,gsm_data[:,f70_gsm],label='gsm')
pylab.scatter(sid_times,data_70*6e10,label='data*1e11',s=1,c='b',edgecolor='b')
pylab.ylim(0,6e3)
pylab.xlim(0,24)
pylab.xlabel('Sidereal Time (Hours)')
pylab.grid()
pylab.savefig(outdir+'gsm_comp_'+sys.argv[2]+'_70_MHz.png',dpi=300)
pylab.clf()

