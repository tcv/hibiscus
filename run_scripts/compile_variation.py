"""
Module to calculate the K_dgsm for a single day of data.
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
import ephem as eph
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('/home/tcv/hibiscus'))
import file_funcs as ff
import eff_funcs as ef
import cal_funcs as cf

#Main directories for the cal and input data
indir = '/lustre/tcv/time_data/'
outdir='/lustre/tcv/time_data/'

date_ind = ['01','03','04','05','06','09','11','12','14']
processed_data=[]
processed_time=[] 
for i in range(0,len(date_ind)):
    filename = indir+'June_'+date_ind[i]+'_Kdgsm_time_series.txt'
    single_data = loadtxt(filename)
    timename = indir+'June_'+date_ind[i]+'_time_indices.txt'
    single_time = loadtxt(timename)
    processed_data.append(single_data)
    processed_time.append(single_time)

#fig = pylab.figure(figsize=(8,6),dpi=300)
#pylab.rc("font",size=6)

pylab.subplot(3,3,1)
pylab.rc("font",size=8)
pylab.scatter(processed_time[0],processed_data[0],c='g',edgecolor='g',s=5)
pylab.xlim(0,24) 
pylab.ylim(0,6000) 
pylab.grid()
pylab.ylabel('Temperature (Kelvin)') 
pylab.legend(['June 1st'],loc=2)

pylab.subplot(3,3,2)
pylab.scatter(processed_time[1],processed_data[1],c='g',edgecolor='g',s=5)
pylab.xlim(0,24) 
pylab.ylim(0,6000) 
pylab.grid()
pylab.legend(['June 3rd'],loc=2)

pylab.subplot(3,3,3)
pylab.scatter(processed_time[2],processed_data[2],c='g',edgecolor='g',s=5)
pylab.xlim(0,24) 
pylab.ylim(0,6000) 
pylab.grid()
pylab.legend(['June 4th'],loc=2)

pylab.subplot(3,3,4) 
pylab.scatter(processed_time[3],processed_data[3],c='g',edgecolor='g',s=5)
pylab.xlim(0,24) 
pylab.ylim(0,6000) 
pylab.grid()
pylab.ylabel('Temperature (Kelvin)')  
pylab.legend(['June 5th'],loc=2)

pylab.subplot(3,3,5) 
pylab.scatter(processed_time[4],processed_data[4],c='g',edgecolor='g',s=5)
pylab.xlim(0,24)  
pylab.ylim(0,6000)  
pylab.grid() 
pylab.legend(['June 6th'],loc=2)
 
pylab.subplot(3,3,6) 
pylab.scatter(processed_time[5],processed_data[5],c='g',edgecolor='g',s=5)
pylab.xlim(0,24)  
pylab.ylim(0,6000)  
pylab.grid() 
pylab.legend(['June 9th'],loc=2)

pylab.subplot(3,3,7)  
pylab.scatter(processed_time[6],processed_data[6],c='g',edgecolor='g',s=5)
pylab.xlim(0,24)  
pylab.ylim(0,6000)  
pylab.grid() 
pylab.ylabel('Temperature (Kelvin)')
pylab.xlabel('Local Sidereal Time (Hours)')  
pylab.legend(['June 11th'],loc=2)
 
pylab.subplot(3,3,8) 
pylab.scatter(processed_time[7],processed_data[7],c='g',edgecolor='g',s=5)
pylab.xlim(0,24)   
pylab.ylim(0,6000)   
pylab.grid()  
pylab.xlabel('Local Sidereal Time (Hours)')  
pylab.legend(['June 12th'],loc=2)
  
pylab.subplot(3,3,9)  
pylab.scatter(processed_time[8],processed_data[8],c='g',edgecolor='g',s=5)
pylab.xlim(0,24)   
pylab.ylim(0,6000)   
pylab.grid()
pylab.xlabel('Local Sidereal Time (Hours)') 
pylab.legend(['June 14th'],loc=2)

pylab.savefig(outdir+'Combined_Kdgsm_time_series',dpi=300) 
pylab.clf() 


