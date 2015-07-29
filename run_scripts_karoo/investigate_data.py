"""
Module for testing out new code to improve flagging and foreground subtraction.
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
import numpy.polynomial.polynomial as poly


indir = sys.argv[1]
outdir = sys.argv[2]

data = numpy.load(indir+'gsm_cal_data_Apr_03_70MHz_ant.npy')
times = numpy.load(indir+'gsm_cal_times_Apr_03_70MHz_ant.npy')
freqs = arange(40.,130.,90./len(data[0]))
f110 = where(freqs<=110.)[0][-1]

mask = zeros((len(data),len(data[0])))
test = where(data==0.0)

mask[test] = 1.0

print 'Percent of Data Masked is: ',100*sum(mask)/(len(mask)*len(mask[0]))

iterate=1
for int in range(0,2):
    mask1 = zeros((len(mask),len(mask[0])))
    mask1[where(mask==1.)]=1.

    bad_ind = zeros(len(freqs))
    for i in range(200,f110):
        if sum(mask[:,i])>0.35*len(mask[:,i]):
             bad_ind[i] = 1.0
#            for f in range(-2,2):
#                bad_ind[i+f]=1.0
    for i in range(0,len(freqs)):
        if bad_ind[i] == 1.:
            mask[:,i] = 1.0
    for i in range(0,len(times)):
        if sum(mask[i])>0.2*len(mask[i]):
            mask[i] = 1.0


print 'Percent of Data Masked after part 1 is: ',100*sum(mask)/(len(mask)*len(mask[0]))

for int in range(0,5):
    
    mask2 = zeros((len(mask),len(mask[0])))
    mask2[where(mask==1.)]=1.
    mean_data = zeros(len(freqs))
    for i in range(0,len(freqs)):
        if len(mask[:,i])==sum(mask[:,i]):
            mean_data[i] = 1.
        else:
            mean = ma.mean(ma.compressed(ma.array(data[:,i],mask=mask[:,i])))
            if mean==0.0:
                mean_data[i] = 1.
            else: 
                mean_data[i] = mean

    mean_mask = zeros(len(mean_data))
    bad_mean = where(mean_data==1.)
    mean_mask[bad_mean] = 1.
    spike_mask = ff.spike_flag(mean_data,mean_mask,freqs,0.5)
    for f in range(0,len(freqs)):
        if spike_mask[f]==1.:
            mask[:,f] = 1.
    f88 = where(freqs<=88.)[0][-1]
    spike_mask_lim = ff.spike_flag(mean_data[0:f88],mean_mask[0:f88],freqs[0:f88],0.2)
    for f in range(0,f88):
        if spike_mask_lim[f]==1.:
            mask[:,f]=1.

    tmean_comp = ma.compressed(ma.array(mean_data,mask=spike_mask))
    freq_comp = ma.compressed(ma.array(freqs,mask=spike_mask))

    print 'Percent of Data Masked after part 1 of iteration ',iterate,' is: ',100*sum(mask)/(len(mask)*len(mask[0]))

    mask3 = zeros((len(mask),len(mask[0])))
    mask3[where(mask==1.)]=1.
    mean_freq = zeros(len(times))
    for i in range(0,len(times)):
        if len(mask[i])==sum(mask[i]):
            mean_freq[i] = 1. 
        else:
            mean = ma.mean(ma.compressed(ma.array(data[i],mask=mask[i])))
            if mean==0.0:
                mean_freq[i] = 1.
            else:
                mean_freq[i] = mean
    
    freq_mask = zeros(len(mean_freq))
    bad_freq = where(mean_freq==1.)
    freq_mask[bad_freq] = 1. 
    spike_freq = ff.timeflag(mean_freq,freq_mask,times,2.,32)
 
    for i in range(0,len(times)):
        if spike_freq[i]==1.:
            mask[i] = 1.

    fmean_comp = ma.compressed(ma.array(mean_freq,mask=spike_freq))
    time_comp = ma.compressed(ma.array(times,mask=spike_freq))

    print 'Percent of Data Masked after part 2 of iteration ',iterate,' is: ',100*sum(mask)/(len(mask)*len(mask[0]))


    pylab.subplot(221)
    pylab.imshow(data,vmin=0,vmax=1.e4,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
    pylab.colorbar()
    pylab.subplot(222)
    pylab.imshow(mask-mask1,vmin=-1,vmax=1.,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
#    pylab.subplot(323)
#    pylab.imshow(mask3-mask2,vmin=-1,vmax=1,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
#    pylab.subplot(324)
#    pylab.imshow(mask-mask3,vmin=-1,vmax=1,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
    pylab.subplot(223)
    pylab.scatter(time_comp,fmean_comp,s=1)
    pylab.subplot(224)
    pylab.scatter(freq_comp,tmean_comp,s=1)
    
    pylab.savefig(outdir+'masking_test_iteration_'+str(iterate)+'.png',dpi=300)
    pylab.clf()

    iterate+=1

for int in range(0,25):
    bad_ind = zeros(len(freqs))
    for i in range(200,f88):
        if sum(mask[:,i])>len(mask[:,i]):
             bad_ind[i] = 1.0
             for f in range(-200,200):
                 bad_ind[i+f]=1.0
    for i in range(0,len(freqs)):
        if bad_ind[i] == 1.:
            mask[:,i] = 1.0
    for i in range(0,len(times)):
        if sum(mask[i])>len(mask[i]):
            mask[i] = 1.0

    mean_data = zeros(len(freqs))
    for i in range(0,len(freqs)):
        if len(mask[:,i])==sum(mask[:,i]):
            mean_data[i] = 1.
        else: 
            mean = ma.mean(ma.compressed(ma.array(data[:,i],mask=mask[:,i])))
            if mean==0.0:
                mean_data[i] = 1.
            else:  
                mean_data[i] = mean

    mean_mask = zeros(len(mean_data))
    bad_mean = where(mean_data==1.)
    mean_mask[bad_mean] = 1.
    spike_mask = ff.spike_flag(mean_data,mean_mask,freqs,1.)
    for f in range(0,len(freqs)):
        if spike_mask[f]==1.: 
            mask[:,f] = 1. 
    tmean_comp = ma.compressed(ma.array(mean_data,mask=spike_mask))
    freq_comp = ma.compressed(ma.array(freqs,mask=spike_mask))
    
print 'Percent of Data Masked after part 3 is: ',100*sum(mask)/(len(mask)*len(mask[0]))


pylab.scatter(freq_comp,tmean_comp,s=1)
pylab.grid()
pylab.xlim(50,110)
pylab.xlabel('Frequency (MHz)')
pylab.ylim(0,6e3)
pylab.ylabel('Temperature (Kelvin)')
pylab.savefig(outdir+'masking_test_final.png',dpi=300)
pylab.clf()

new_data = data
masked = where(mask==1.)
new_data[masked]=0.0

numpy.save(indir+'gsm_cal_data_masked_Apr_03_70MHz_ant_shortsub.npy',new_data)
