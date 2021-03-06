"""
Module for testing out new code to improve flagging.
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
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
direct = sys.argv[3]


data = np.load(indir+'gsm_cal_data_Apr_'+direct+'MHz_ant.npy')
times = np.load(indir+'gsm_cal_times_Apr_'+direct+'MHz_ant.npy')
freqs = np.arange(40.,130.,90./len(data[0]))
f110 = np.where(freqs<=110.)[0][-1]
f90 = np.where(freqs<=90.)[0][-1]
mask = np.zeros((len(data),len(data[0])))
test = np.where(data==0.0)

mask[test] = 1.0

if direct=='06_100':
    for t in range(0,len(times)):
        if times[t]>5:
            if times[t]<13:
                mask[t] = np.ones(len(freqs))
                data[t] = np.zeros(len(freqs))
print 'Percent of Data Masked is: ',100*np.sum(mask)/(len(mask)*len(mask[0]))

iterate=1
for int in range(0,2):
    mask1 = np.zeros((len(mask),len(mask[0])))
    mask1[np.where(mask==1.)]=1.

    bad_ind = np.zeros(len(freqs))
    for i in range(f90,f110):
        if sum(mask[:,i])>0.5*len(mask[:,i]):
             bad_ind[i] = 1.0
#            for f in range(-2,2):
#                bad_ind[i+f]=1.0
    for i in range(0,len(freqs)):
        if bad_ind[i] == 1.:
            mask[:,i] = 1.0
    for i in range(0,len(times)):
        if sum(mask[i])>0.5*len(mask[i]):
            mask[i] = 1.0


print 'Percent of Data Masked after part 1 is: ',100*np.sum(mask)/(len(mask)*len(mask[0]))

for int in range(0,5):
    
    mask2 = np.zeros((len(mask),len(mask[0])))
    mask2[np.where(mask==1.)]=1.
    mean_data = np.zeros(len(freqs))
    for i in range(0,len(freqs)):
        if len(mask[:,i])==sum(mask[:,i]):
            mean_data[i] = 1.
        else:
            mean = ma.mean(ma.compressed(ma.array(data[:,i],mask=mask[:,i])))
            if mean==0.0:
                mean_data[i] = 1.
            else: 
                mean_data[i] = mean

    mean_mask = np.zeros(len(mean_data))
    bad_mean = np.where(mean_data==1.)
    mean_mask[bad_mean] = 1.
    spike_mask = ff.spike_flag(mean_data,mean_mask,freqs,1.)
    mean_mask = spike_mask
    for f in range(0,len(freqs)):
        if spike_mask[f]==1.:
            mask[:,f] = 1.
    if sys.argv[3].split('_')[-1]=='100':
        finit = np.where(freqs<=50.)[0][-1]
        ffinal = np.where(freqs<=110.)[0][-1]
    elif sys.argv[3].split('_')[-1]=='70':
        finit = np.where(freqs<=50.)[0][-1]
        ffinal = np.where(freqs<=110.)[0][-1]
    f88 = np.where(freqs<=88.)[0][-1]

    spike_mask_lim = ff.spike_flag(mean_data[finit:ffinal],mean_mask[finit:ffinal],freqs[finit:ffinal],1.)
    for f in range(0,ffinal-finit):
        if spike_mask_lim[f]==1.:
            mask[:,f+finit]=1.

    tmean_comp = ma.compressed(ma.array(mean_data,mask=mean_mask))
    freq_comp = ma.compressed(ma.array(freqs,mask=mean_mask))

    print 'Percent of Data Masked after part 1 of iteration ',iterate,' is: ',100*np.sum(mask)/(len(mask)*len(mask[0]))

    mask3 = np.zeros((len(mask),len(mask[0])))
    mask3[np.where(mask==1.)]=1.
    mean_freq = np.zeros(len(times))
    for i in range(0,len(times)):
        if len(mask[i])==sum(mask[i]):
            mean_freq[i] = 1. 
        else:
            mean = ma.mean(ma.compressed(ma.array(data[i],mask=mask[i])))
            if mean==0.0:
                mean_freq[i] = 1.
            else:
                mean_freq[i] = mean
    
    freq_mask = np.zeros(len(mean_freq))
    bad_freq = np.where(mean_freq==1.)
    freq_mask[bad_freq] = 1. 
    spike_freq = ff.timeflag(mean_freq,freq_mask,times,2.,32)
 
    for i in range(0,len(times)):
        if spike_freq[i]==1.:
            mask[i] = 1.

    fmean_comp = ma.compressed(ma.array(mean_freq,mask=spike_freq))
    time_comp = ma.compressed(ma.array(times,mask=spike_freq))

    print 'Percent of Data Masked after part 2 of iteration ',iterate,' is: ',100*np.sum(mask)/(len(mask)*len(mask[0]))

    pylab.rc('font',size=6)
    pylab.subplot(221)
    pylab.imshow(data,vmin=0,vmax=6.e3,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
    pylab.colorbar()
    pylab.xlabel('Sidereal Time (Hours)')
    pylab.ylabel('Frequency (MHz)')
    pylab.title('Calibrated Data')
    pylab.subplot(222)
    pylab.imshow(mask-mask1,vmin=-.1,vmax=1.1,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
    pylab.title('Current Mask-Initial Mask')
    pylab.xlabel('Sidereal Time (Hours)')
    pylab.ylabel('Frequency (MHz)') 
    pylab.subplot(223)
    pylab.scatter(time_comp,fmean_comp/1000.,s=1)
    pylab.ylim(0,6)
    pylab.xlim(0,24)
    pylab.ylabel('Temperature (1000 K)')
    pylab.xlabel('Sidereal Time (Hours)')
    pylab.subplot(224)
    pylab.scatter(freq_comp,tmean_comp/1000.,s=1)
    pylab.xlim(40,130)
    pylab.ylim(0,6.)
    pylab.xlabel('Frequency (MHz)')    
    pylab.savefig(outdir+'masking_test_iteration_'+str(iterate)+'.png',dpi=300)
    pylab.clf()

    iterate+=1

for int in range(0,25):
    bad_ind = np.zeros(len(freqs))
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

    mean_data = np.zeros(len(freqs))
    for i in range(0,len(freqs)):
        if len(mask[:,i])==sum(mask[:,i]):
            mean_data[i] = 1.
        else: 
            mean = ma.mean(ma.compressed(ma.array(data[:,i],mask=mask[:,i])))
            if mean==0.0:
                mean_data[i] = 1.
            else:  
                mean_data[i] = mean

    mean_mask = np.zeros(len(mean_data))
    bad_mean = np.where(mean_data==1.)
    mean_mask[bad_mean] = 1.
    spike_mask = ff.spike_flag(mean_data,mean_mask,freqs,1.)
    for f in range(0,len(freqs)):
        if spike_mask[f]==1.: 
            mask[:,f] = 1. 
    tmean_comp = ma.compressed(ma.array(mean_data,mask=spike_mask))
    freq_comp = ma.compressed(ma.array(freqs,mask=spike_mask))
    
print 'Percent of Data Masked after part 3 is: ',100*np.sum(mask)/(len(mask)*len(mask[0]))

pylab.rc('font',size=8)
f,axarr = plt.subplots(1,2)
axarr[0].scatter(freq_comp,tmean_comp,s=1)
axarr[0].grid()
axarr[0].set_xlim(50,110)
axarr[0].set_xlabel('Frequency (MHz)')
axarr[0].set_ylim(1e3,6e3)
axarr[0].set_ylabel('Temperature (Kelvin)')
axarr[0].set_title('Mean Data')
axarr[1].imshow(mask,vmin=0,vmax=1.,aspect=90./24.,extent=(freqs[0],freqs[-1],24.,0))
axarr[1].set_title('Full Mask')
pylab.savefig(outdir+'masking_test_final.png',dpi=300)
pylab.clf()

new_data = data
masked = np.where(mask==1.)
new_data[masked]=0.0

np.save(indir+'gsm_cal_data_masked_Apr_'+direct+'MHz_ant.npy',new_data)
