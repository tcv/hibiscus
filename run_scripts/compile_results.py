"""
Module to remove foreground signal using polynomial fitting
for a single day of data. 
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
import numpy.polynomial.polynomial as poly
sys.path.append(os.path.abspath('/home/tcv/hibiscus'))
import file_funcs as ff
import eff_funcs as ef
import cal_funcs as cf

#Main directories for the input and output
indir = '/lustre/tcv/mean_cal_data/'
outdir = '/lustre/tcv/mean_cal_data/'
directories = os.listdir(indir)

# Make data array
Kt_data = []
Kt_mask = []
Kt_fit = []
Kdgsm_data = []
Kdgsm_mask = []
Kdgsm_fit = []
Kgsm_data = []
Kgsm_mask = []
Kgsm_fit = []
dates = arange(1,15)
#limit dates to remove bad days (02,05,08,10)
dates = [1,3,4,6,7,9,11,12,13,14]

for date in dates:
    if date<10:
        str_date = '0'+str(date)
    else:
        str_date = str(date)
    sKt_data = loadtxt(indir+'June_'+str_date+'_Kt_mean.txt')
    Kt_data.append(sKt_data)
    sKt_mask = loadtxt(indir+'June_'+str_date+'_Kt_mean_mask.txt')
    Kt_mask.append(sKt_mask)
    sKt_fit = loadtxt(indir+'June_'+str_date+'_Kt_fit2.txt')
    Kt_fit.append(sKt_fit)
    sKdgsm_data = loadtxt(indir+'June_'+str_date+'_Kdgsm_mean.txt')
    Kdgsm_data.append(sKdgsm_data)
    sKdgsm_mask = loadtxt(indir+'June_'+str_date+'_Kdgsm_mean_mask.txt')
    Kdgsm_mask.append(sKdgsm_mask)
    sKdgsm_fit = loadtxt(indir+'June_'+str_date+'_Kdgsm_fit2.txt')
    Kdgsm_fit.append(sKdgsm_fit)
    sKgsm_data = loadtxt(indir+'June_'+str_date+'_Kgsm_mean.txt')
    Kgsm_data.append(sKgsm_data)
    sKgsm_mask = loadtxt(indir+'June_'+str_date+'_Kgsm_mean_mask.txt')
    Kgsm_mask.append(sKgsm_mask)
    sKgsm_fit = loadtxt(indir+'June_'+str_date+'_Kgsm_fit2.txt')
    Kgsm_fit.append(sKgsm_fit)

Kt_data = array(Kt_data)
Kt_mask = array(Kt_mask)
Kt_fit = array(Kt_fit)
Kdgsm_data = array(Kdgsm_data) 
Kdgsm_mask = array(Kdgsm_mask) 
Kdgsm_fit = array(Kdgsm_fit)
Kgsm_data = array(Kgsm_data) 
Kgsm_mask = array(Kgsm_mask) 
Kgsm_fit = array(Kgsm_fit)

width = 90.0/len(sKt_data)
freq = arange(40.+width/2,130.,width)

#Before rebinning in frequency
Kt_resid = ma.array(Kt_data-Kt_fit,mask=Kt_mask)
Kdgsm_resid =ma.array( Kdgsm_data-Kdgsm_fit,mask=Kdgsm_mask)
Kgsm_resid = ma.array(Kgsm_data-Kgsm_fit,mask=Kgsm_mask)

#Individual daily residuals plotting
plot_c = ['b','g','r','c','m','y','k']
for i in range(0,len(Kt_resid)/2):
    pylab.scatter(freq,Kt_resid[i],c=plot_c[i],s=10,label=str(dates[i]))
pylab.xlim(50.,90.)
pylab.ylim(-100,100)
pylab.grid()
pylab.legend()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Daily Residuals from log fits - Kt cal')
pylab.savefig(outdir+'June_Kt_residuals_pt1',dpi=300)
pylab.clf()
for i in range(len(Kt_resid)/2,len(Kt_resid)):
    pylab.scatter(freq,Kt_resid[i],c=plot_c[i-7],s=10,label=str(dates[i]))
pylab.xlim(50.,90.) 
pylab.ylim(-100,100) 
pylab.grid() 
pylab.legend() 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Daily Residuals from log fits - Kt cal')
pylab.savefig(outdir+'June_Kt_residuals_pt2',dpi=300)
pylab.clf() 

for i in range(0,len(Kdgsm_resid)/2):
    pylab.scatter(freq,Kdgsm_resid[i],c=plot_c[i],s=10,label=str(dates[i]))
pylab.xlim(50.,90.) 
pylab.ylim(-100,100) 
pylab.grid() 
pylab.legend() 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Residuals (Kelvin)') 
pylab.title('Daily Residuals from log fits - Kdgsm cal')
pylab.savefig(outdir+'June_Kdgsm_residuals_pt1',dpi=300)
pylab.clf() 
for i in range(len(Kdgsm_resid)/2,len(Kdgsm_resid)):
    pylab.scatter(freq,Kdgsm_resid[i],c=plot_c[i-7],s=10,label=str(dates[i]))
pylab.xlim(50.,90.)  
pylab.ylim(-100,100)  
pylab.grid()  
pylab.legend()  
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Residuals (Kelvin)') 
pylab.title('Daily Residuals from log fits - Kdgsm cal')
pylab.savefig(outdir+'June_Kdgsm_residuals_pt2',dpi=300)
pylab.clf()  

for i in range(0,len(Kgsm_resid)/2):
    pylab.scatter(freq,Kgsm_resid[i],c=plot_c[i],s=10,label=str(dates[i]))
pylab.xlim(50.,90.) 
pylab.ylim(-10,10) 
pylab.grid() 
pylab.legend() 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Residuals (Kelvin)') 
pylab.title('Daily Residuals from log fits - Kgsm cal')
pylab.savefig(outdir+'June_Kgsm_residuals_pt1',dpi=300)
pylab.clf() 
for i in range(len(Kgsm_resid)/2,len(Kgsm_resid)):
    pylab.scatter(freq,Kgsm_resid[i],c=plot_c[i-7],s=10,label=str(dates[i]))
pylab.xlim(50.,90.)  
pylab.ylim(-10,10)  
pylab.grid()  
pylab.legend()  
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Residuals (Kelvin)') 
pylab.title('Daily Residuals from log fits - Kgsm cal')
pylab.savefig(outdir+'June_Kgsm_residuals_pt2',dpi=300)
pylab.clf()  


Kt_rmean = ma.mean(Kt_resid,axis=0)
Kdgsm_rmean = ma.mean(Kdgsm_resid,axis=0)
Kgsm_rmean = ma.mean(Kgsm_resid,axis=0)
Kt_rstd = ma.std(Kt_resid,axis=0)
Kdgsm_rstd = ma.std(Kdgsm_resid,axis=0)
Kgsm_rstd = ma.std(Kgsm_resid,axis=0)

pylab.scatter(freq,Kt_rmean,c='b',edgecolor='b',s=5,label='Kt mean')
pylab.errorbar(freq,Kt_rmean,Kt_rstd,c='b',fmt=None,label='Kt std')
pylab.scatter(freq,Kdgsm_rmean,c='g',edgecolor='g',s=5,label='Kdgsm mean')
pylab.errorbar(freq,Kdgsm_rmean,Kdgsm_rstd,c='b',fmt=None,label='Kdgsm std')
pylab.scatter(freq,Kgsm_rmean,c='r',edgecolor='r',s=5,label='Kgsm mean')
pylab.errorbar(freq,Kgsm_rmean,Kgsm_rstd,c='r',fmt=None,label='Kgsm std')
pylab.xlim(50.,90.)
pylab.ylim(-100.,100)
pylab.grid()
pylab.legend()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits - No modifications')
pylab.savefig(outdir+'June_original_residuals',dpi=300)
pylab.clf()

Kt_rbdata = []
Kt_rbmask = []
Kdgsm_rbdata = []
Kdgsm_rbmask = []
Kgsm_rbdata = []
Kgsm_rbmask = []

for i in range(0,len(dates)):
    Kt_data_rb,Kt_mask_rb,freq_rb = ff.rebin(Kt_data[i]-Kt_fit[i],Kt_mask[i],freq,10)
    Kt_rbdata.append(Kt_data_rb)
    Kt_rbmask.append(Kt_mask_rb)
    Kdgsm_data_rb,Kdgsm_mask_rb,freq_rb = ff.rebin(Kdgsm_data[i]-Kdgsm_fit[i],Kdgsm_mask[i],freq,10)
    Kdgsm_rbdata.append(Kdgsm_data_rb)
    Kdgsm_rbmask.append(Kdgsm_mask_rb)
    Kgsm_data_rb,Kgsm_mask_rb,freq_rb = ff.rebin(Kgsm_data[i]-Kgsm_fit[i],Kgsm_mask[i],freq,10)
    Kgsm_rbdata.append(Kgsm_data_rb)
    Kgsm_rbmask.append(Kgsm_mask_rb)

Kt_rbdata = ma.array(Kt_rbdata,mask=Kt_rbmask)
Kdgsm_rbdata = ma.array(Kdgsm_rbdata,mask=Kdgsm_rbmask)
Kgsm_rbdata = ma.array(Kgsm_rbdata,mask=Kgsm_rbmask)

Kt_rbmean = ma.mean(Kt_rbdata,axis=0)
Kdgsm_rbmean = ma.mean(Kdgsm_rbdata,axis=0)
Kgsm_rbmean = ma.mean(Kgsm_rbdata,axis=0)
Kt_rbstd = ma.std(Kt_rbdata,axis=0)
Kdgsm_rbstd = ma.std(Kdgsm_rbdata,axis=0)
Kgsm_rbstd = ma.std(Kgsm_rbdata,axis=0)

pylab.scatter(freq_rb,Kt_rbmean,c='b',edgecolor='b',s=5,label='Kt mean')
pylab.errorbar(freq_rb,Kt_rbmean,Kt_rbstd,c='b',fmt=None,label='Kt std') 
pylab.scatter(freq_rb,Kdgsm_rbmean,c='g',edgecolor='g',s=5,label='Kdgsm mean')
pylab.errorbar(freq_rb,Kdgsm_rbmean,Kdgsm_rbstd,c='b',fmt=None,label='Kdgsm std')
pylab.scatter(freq_rb,Kgsm_rbmean,c='r',edgecolor='r',s=5,label='Kgsm mean')
pylab.errorbar(freq_rb,Kgsm_rbmean,Kgsm_rbstd,c='r',fmt=None,label='Kgsm std')
pylab.xlim(50.,90.) 
pylab.ylim(-50.,50) 
pylab.grid() 
pylab.legend()
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits - rebinned')
pylab.savefig(outdir+'June_rebinned_residuals',dpi=300)
pylab.clf()




