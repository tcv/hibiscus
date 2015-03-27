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

#date_ind = ['01','03','04','05','06','09','11','12','14']
date_ind = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
times = arange(0,24,24./289)

psingle_Kdgsm_data=[]
psingle_Kt_data=[]
p_time=[] 
pfull_Kdgsm_data= zeros((len(date_ind),len(times),205))
pfull_Kt_data = zeros((len(date_ind),len(times),205))
p_mask_data_Kt = ones((len(date_ind),len(times),205))
p_mask_data_Kdgsm = ones((len(date_ind),len(times),205))


#for i in range(0,len(date_ind)):
for i in range(0,len(date_ind)):
    filename = indir+'June_'+date_ind[i]+'_Kdgsm_time_series.txt'
    single_data = loadtxt(filename)
    psingle_Kdgsm_data.append(single_data)
    filename = indir+'June_'+date_ind[i]+'_Kt_time_series.txt'
    single_data = loadtxt(filename)
    psingle_Kt_data.append(single_data)
    timename = indir+'June_'+date_ind[i]+'_time_indices.txt'
    single_time = loadtxt(timename)
    p_time.append(single_time)
    filename = indir+'June_'+date_ind[i]+'_Kdgsm_full_time_series.txt'
    single_data = loadtxt(filename)
    print shape(single_data)
    print where(times==single_time[0])[0][0]
    for t in range(0,len(single_time)):
        time_ind = where(times==single_time[t])[0][0]
        pfull_Kdgsm_data[i,time_ind] = single_data[t]

    filename = indir+'June_'+date_ind[i]+'_Kt_full_time_series.txt'
    single_data = loadtxt(filename)
    for t in range(0,len(single_time)):
        time_ind = where(times==single_time[t])[0][0]
#        print time_ind
        pfull_Kt_data[i,time_ind] = single_data[t]

#    pfull_Kt_data.append(single_data) 
    filename = indir+'June_'+date_ind[i]+'_mask_full_time_series.txt'
    single_data = loadtxt(filename)
    for t in range(0,len(single_time)):
        time_ind = where(times==single_time[t])[0][0]
        p_mask_data_Kt[i,time_ind] = single_data[t]
        p_mask_data_Kdgsm[i,time_ind] = single_data[t]

#    p_mask_data.append(single_data) 
    
freqs = arange(50,100.,50./205)
print len(freqs)

#Avoid bad freq ranges set by looking at avg cal plots
for i in range(0,len(date_ind)):
    if date_ind[i]=='01':
        tbad_min = where(times<=12.)[0][-1]
        tbad_max = where(times<=14.5)[0][-1]
        for f in range(0,len(freqs)):
            for t in range(0,len(times)):
                p_mask_data_Kt[i,t,f] = 1.0
#            for t in range(tbad_min,tbad_max):
                p_mask_data_Kdgsm[i,t,f] = 1.0


    if date_ind[i]=='02':
        fbad_min = where(freqs<=55.)[0][-1]
        fbad_max = where(freqs<=60.)[0][-1]
        for t in range(0,len(times)):
            for f in range(0,len(freqs)):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
    
    if date_ind[i]=='03':
        tbad_min = where(times<=1.4)[0][-1]
        tbad_max = where(times<=9.4)[0][-1]
        for t in range(tbad_min,tbad_max):
            for f in range(0,len(freqs)):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='04':
        fbad_min = where(freqs<=78.)[0][-1]
        fbad_max = where(freqs<=84.)[0][-1]
        tbad_min = where(times<=3.)[0][-1]
        tbad_max = where(times<=6.0)[0][-1]
        for t in range(tbad_min,tbad_max):
            for f in range(0,len(freqs)):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='05':
        for t in range(0,len(times)):
            for f in range(0,len(freqs)):
                p_mask_data_Kdgsm[i,t,f] = 1.0
   
    if date_ind[i]=='06':
        fbad_min= where(freqs<=59.5)[0][-1]
        fbad_max = where(freqs<=71.)[0][-1]
        fbad_min2 = where(freqs<=76.)[0][-1]
        fbad_max2 = where(freqs<=84.)[0][-1]
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
            for f in range(fbad_min2,fbad_max2):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='07':
        fbad_min = where(freqs<=59.)[0][-1]
        fbad_max = where(freqs<=71.)[0][-1]
        fbad_min2 = where(freqs<=75.)[0][-1]
        fbad_max2 = where(freqs<=83.)[0][-1]
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
            for f in range(fbad_min2,fbad_max2):
                p_mask_data_Kt[i,t,f] = 1.0
            for f in range(0,len(freqs)):
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='08':
        fbad_min = where(freqs<=55.)[0][-1]
        fbad_max = where(freqs<=64.)[0][-1]
        fbad_min2 = where(freqs<=76.)[0][-1]
        fbad_max2 = where(freqs<=81.)[0][-1]
        tbad_min = where(times<=2.)[0][-1]
        tbad_max = where(times<=5.)[0][-1]
        for t in range(tbad_min,tbad_max):
            for f in range(0,len(freqs)):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
            for f in range(fbad_min2,fbad_max2):
                p_mask_data_Kt[i,t,f] = 1.0
            for f in range(0,len(freqs)):
                p_mask_data_Kdgsm[i,t,f] = 1.0
                
    if date_ind[i]=='09':
        fbad_min = where(freqs<=55.)[0][-1]
        fbad_max = where(freqs<=65.)[0][-1]
        fbad_min2 = where(freqs<=77.)[0][-1]
        fbad_max2 = where(freqs<=81.)[0][-1]
        tbad_min = where(times<=19.5)[0][-1]
        tbad_max = where(times<=24.)[0][-1]
        for t in range(tbad_min,tbad_max):
            for f in range(0,len(freqs)):
                p_mask_data_Kt[i,t,f] = 1.0
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
            for f in range(fbad_min2,fbad_max2):
                p_mask_data_Kt[i,t,f] = 1.0
            for f in range(0,len(freqs)):
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='10':
        for f in range(0,len(freqs)):
            for t in range(0,len(times)):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='11':
        fbad_min = where(freqs<=65.)[0][-1]
        fbad_max = where(freqs<=72.)[0][-1]
        fbad_min2 = where(freqs<=80.)[0][-1]
        fbad_max2 = where(freqs<=88.)[0][-1]
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
            for f in range(fbad_min2,fbad_max2):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='12':
        fbad_min = where(freqs<=55.)[0][-1]
        fbad_max = where(freqs<=62.)[0][-1]
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0

    if date_ind[i]=='13':
        fbad_min = where(freqs<=55.)[0][-1]
        fbad_max = where(freqs<=62.)[0][-1]
        for f in range(fbad_min,fbad_max):
            for t in range(0,len(times)):
                p_mask_data_Kt[i,t,f] = 1.0
        for f in range(0,len(freqs)):
            for t in range(0,len(times)): 
                p_mask_data_Kdgsm[i,t,f] = 1.0
                
    if date_ind[i]=='14':
        fbad_min = where(freqs<=78.)[0][-1]
        fbad_max = where(freqs<=84.)[0][-1]
        tbad_min = where(times<=9.)[0][-1]
        tbad_max = where(times<=12.)[0][-1]
        for t in range(tbad_min,tbad_max):
            for f in range(0,len(freqs)):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
        for t in range(0,len(times)):
            for f in range(fbad_min,fbad_max):
                p_mask_data_Kt[i,t,f] = 1.0
                p_mask_data_Kdgsm[i,t,f] = 1.0
        
print 'Percentage of Kt data masked',100*sum(p_mask_data_Kt)/(len(p_mask_data_Kt)*len(p_mask_data_Kt[0])*len(p_mask_data_Kt[0,0]))
print 'Percentage of Kdgsm data masked',100*sum(p_mask_data_Kdgsm)/(len(p_mask_data_Kt)*len(p_mask_data_Kt[0])*len(p_mask_data_Kt[0,0]))


binscale = 7
prebin_Kdgsm_data = zeros((len(date_ind),len(times),int(205/binscale)))
prebin_Kt_data = zeros((len(date_ind),len(times),int(205/binscale)))
prebin_mask_Kt = zeros((len(date_ind),len(times),int(205/binscale))) 
prebin_mask_Kdgsm = zeros((len(date_ind),len(times),int(205/binscale)))
print int(205/binscale)
for t in range(0,len(times)):
    for ind in range(0,len(date_ind)):
        new_data,new_mask,new_freq = ff.rebin(pfull_Kdgsm_data[ind,t],p_mask_data_Kdgsm[ind,t],freqs,binscale)
#        #print len(new_data)
        prebin_Kdgsm_data[ind,t] = new_data
        prebin_mask_Kdgsm[ind,t] = new_mask
        new_data,new_mask,new_freq = ff.rebin(pfull_Kt_data[ind,t],p_mask_data_Kt[ind,t],freqs,binscale)
        prebin_Kt_data[ind,t] = new_data
        prebin_mask_Kt[ind,t] = new_mask

Kt_fit_full = zeros((len(date_ind),len(times),len(new_freq)))
Kt_params_full = ones((len(date_ind),len(times),3))
Kdgsm_fit_full = zeros((len(date_ind),len(times),len(new_freq)))
Kdgsm_params_full = ones((len(date_ind),len(times),3))
stds = ones(len(new_freq))
for d in range(0,len(date_ind)):
    for t in range(0,len(times)):
        if len(ma.compressed(ma.array(new_freq,mask=prebin_mask_Kt[d,t])))>1.:
            Kt_fit_full[d,t],Kt_params_full[d,t] = cf.poly_fore(prebin_Kt_data[d,t],prebin_mask_Kt[d,t],new_freq,60.,90.,2,stds)
        if len(ma.compressed(ma.array(new_freq,mask=prebin_mask_Kdgsm[d,t])))>1.: 
            Kdgsm_fit_full[d,t],Kdgsm_params_full[d,t]=cf.poly_fore(prebin_Kdgsm_data[d,t],prebin_mask_Kdgsm[d,t],new_freq,60.,90.,2,stds)

Kt_params_mean = ones((len(times),3))
Kdgsm_params_mean = ones((len(times),3))
Kt_params_std = zeros((len(times),3))
Kdgsm_params_std = zeros((len(times),3))

for t in range(0,len(times)):
    Kt_mask_ind = where(Kt_params_full[:,t,0]!=1)[0]
    Kt_params_mean[t,0] = ma.mean(Kt_params_full[Kt_mask_ind,t,0])
    Kt_params_mean[t,1] = ma.mean(Kt_params_full[Kt_mask_ind,t,1])
    Kt_params_mean[t,2] = ma.mean(Kt_params_full[Kt_mask_ind,t,2])
    Kdgsm_mask_ind = where(Kdgsm_params_full[:,t,0]!=1)[0]
    Kdgsm_params_mean[t,0] = ma.mean(Kdgsm_params_full[Kdgsm_mask_ind,t,0])
    Kdgsm_params_mean[t,1] = ma.mean(Kdgsm_params_full[Kdgsm_mask_ind,t,1])
    Kdgsm_params_mean[t,2] = ma.mean(Kdgsm_params_full[Kdgsm_mask_ind,t,2])
    Kt_params_std[t,0] = ma.std(Kt_params_full[Kt_mask_ind,t,0])
    Kt_params_std[t,1] = ma.std(Kt_params_full[Kt_mask_ind,t,1])
    Kt_params_std[t,2] = ma.std(Kt_params_full[Kt_mask_ind,t,2])
    Kdgsm_params_std[t,0] = ma.std(Kdgsm_params_full[Kdgsm_mask_ind,t,0])
    Kdgsm_params_std[t,1] = ma.std(Kdgsm_params_full[Kdgsm_mask_ind,t,1])
    Kdgsm_params_std[t,2] = ma.std(Kdgsm_params_full[Kdgsm_mask_ind,t,2])


pylab.rc("font",size=7)
pylab.subplot(1,3,1)
pylab.scatter(times,Kt_params_mean[:,0],c='b',edgecolor='b',s=5)
pylab.scatter(times,Kdgsm_params_mean[:,0],c='g',edgecolor='g',s=5)
#pylab.errorbar(times,Kt_params_mean[:,0],Kt_params_std[:,0],ecolor='b',fmt=None)
#pylab.errorbar(times,Kdgsm_params_mean[:,0],Kdgsm_params_std[:,0],ecolor='g',fmt=None)
pylab.xlim(0,24)
pylab.xlabel('LST (Hours)')
pylab.ylim(3.1,3.8)
pylab.legend(('Kt','Kdgsm'))
pylab.grid()

pylab.subplot(1,3,2)
pylab.scatter(times,Kt_params_mean[:,1],c='b',edgecolor='b',s=5)
pylab.scatter(times,Kdgsm_params_mean[:,1],c='g',edgecolor='g',s=5)
#pylab.errorbar(times,Kt_params_mean[:,1],Kt_params_std[:,1],ecolor='b',fmt=None)
#pylab.errorbar(times,Kdgsm_params_mean[:,1],Kdgsm_params_std[:,1],ecolor='g',fmt=None) 
pylab.xlim(0,24)
pylab.xlabel('LST (Hours)')
pylab.ylim(-2.7,-2.1)
pylab.legend(('Kt','Kdgsm'))
pylab.grid()

pylab.subplot(1,3,3)
pylab.scatter(times,Kt_params_mean[:,2],c='b',edgecolor='b',s=5)
pylab.scatter(times,Kdgsm_params_mean[:,2],c='g',edgecolor='g',s=5)
#pylab.errorbar(times,Kt_params_mean[:,2],Kt_params_std[:,2],ecolor='b',fmt=None)
#pylab.errorbar(times,Kdgsm_params_mean[:,2],Kdgsm_params_std[:,2],ecolor='g',fmt=None) 
pylab.xlim(0,24)
pylab.xlabel('LST (Hours)')
pylab.ylim(-5,5)
pylab.legend(('Kt','Kdgsm'))
pylab.grid()

pylab.savefig(outdir+'Combined_mean_fore_params',dpi=300)
pylab.clf()    

pylab.rc("font",size=7)
pylab.subplot(2,3,1)
pylab.imshow(Kt_params_full[:,:,0],vmin=3.1,vmax=3.7,aspect=(times[-1]-times[0])/14,extent=(times[0],times[-1],14,1))
pylab.colorbar()
pylab.ylabel('Day in June 2013')
pylab.title('Kt a0')

pylab.subplot(2,3,2)
pylab.imshow(Kt_params_full[:,:,1],vmin=-2.85,vmax=-2.0,aspect=(times[-1]-times[0])/14,extent=(times[0],times[-1],14,1))
pylab.colorbar()
pylab.title('Kt a1')

pylab.subplot(2,3,3) 
pylab.imshow(Kt_params_full[:,:,2],vmin=-5,vmax=5,aspect=(times[-1]-times[0])/14,extent=(times[0],times[-1],14,1))
pylab.colorbar()
pylab.title('Kt a2')

pylab.subplot(2,3,4)
pylab.imshow(Kdgsm_params_full[:,:,0],vmin=3.1,vmax=3.7,aspect=(times[-1]-times[0])/14,extent=(times[0],times[-1],14,1))
pylab.colorbar()
pylab.ylabel('Day in June 2013')
pylab.title('Kdgsm a0')
pylab.xlabel('LST (Hours)')

pylab.subplot(2,3,5) 
pylab.imshow(Kdgsm_params_full[:,:,1],vmin=-2.85,vmax=-2.0,aspect=(times[-1]-times[0])/14,extent=(times[0],times[-1],14,1))
pylab.colorbar()
pylab.title('Kdgsm a1')
pylab.xlabel('LST (Hours)')

pylab.subplot(2,3,6)
pylab.imshow(Kdgsm_params_full[:,:,2],vmin=-5,vmax=5,aspect=(times[-1]-times[0])/14,extent=(times[0],times[-1],14,1))
pylab.colorbar() 
pylab.title('Kdgsm a2')
pylab.xlabel('LST (Hours)')

pylab.savefig(outdir+'Combined_fore_params',dpi=300)
pylab.clf()


p_tmean_Kdgsm_data = zeros((len(times),len(new_freq)))
p_tstd_Kdgsm_data = zeros((len(times),len(new_freq)))
p_tmean_Kt_data = zeros((len(times),len(new_freq)))
p_tstd_Kt_data = zeros((len(times),len(new_freq)))
p_tmask_Kt = ones((len(times),len(new_freq)))
p_tmask_Kdgsm = ones((len(times),len(new_freq)))
p_tdiff_Kdgsm_data = zeros((len(times),len(new_freq)))
p_tdiff_Kt_data = zeros((len(times),len(new_freq)))
p_tdstd_Kdgsm_data = zeros((len(times),len(new_freq)))
p_tdstd_Kt_data = zeros((len(times),len(new_freq)))

for f in range(0,len(new_freq)):
    for t in range(0,len(times)):
        if len(prebin_mask_Kdgsm[:,t,f])>sum(prebin_mask_Kdgsm[:,t,f]):
            single_mean = ma.mean(ma.compressed(ma.array(prebin_Kdgsm_data[:,t,f],mask=prebin_mask_Kdgsm[:,t,f])))
            single_std = ma.std(ma.compressed(ma.array(prebin_Kdgsm_data[:,t,f],mask=prebin_mask_Kdgsm[:,t,f])))
            if isnan(single_mean)==True:
                p_tmask_Kdgsm[t,f]=1.0
                p_tmean_Kdgsm_data[t,f] = 0.
                p_tstd_Kdgsm_data[t,f] = 0.
            else:
                p_tmask_Kdgsm[t,f] = 0.
                p_tmean_Kdgsm_data[t,f] = single_mean
                p_tstd_Kdgsm_data[t,f] = single_std
            single_mean = ma.mean(ma.compressed(ma.array(prebin_Kdgsm_data[:,t,f]-Kdgsm_fit_full[:,t,f],mask=prebin_mask_Kdgsm[:,t,f])))
            single_std = ma.std(ma.compressed(ma.array(prebin_Kdgsm_data[:,t,f]-Kdgsm_fit_full[:,t,f],mask=prebin_mask_Kdgsm[:,t,f])))
            if isnan(single_mean)==True:
                p_tdiff_Kdgsm_data[t,f] = 0.
                p_tdstd_Kdgsm_data[t,f] = 0.
            else:
                p_tdiff_Kdgsm_data[t,f] = single_mean
                p_tdstd_Kdgsm_data[t,f] = single_std

        if len(prebin_mask_Kt[:,t,f])>sum(prebin_mask_Kt[:,t,f]): 
            single_mean = ma.mean(ma.compressed(ma.array(prebin_Kt_data[:,t,f],mask=prebin_mask_Kt[:,t,f])))
            single_std = ma.std(ma.compressed(ma.array(prebin_Kt_data[:,t,f],mask=prebin_mask_Kt[:,t,f])))
            if isnan(single_mean)==True:
                p_tmask_Kt[t,f] = 1.0
                p_tmean_Kt_data[t,f] = 0.
                p_tstd_Kt_data[t,f] = 0.
            else:
                p_tmask_Kt[t,f] = 0.
                p_tmean_Kt_data[t,f] = single_mean
                p_tstd_Kt_data[t,f] = single_std
            single_mean = ma.mean(ma.compressed(ma.array(prebin_Kt_data[:,t,f]-Kt_fit_full[:,t,f],mask=prebin_mask_Kt[:,t,f])))
            single_std = ma.std(ma.compressed(ma.array(prebin_Kt_data[:,t,f]-Kt_fit_full[:,t,f],mask=prebin_mask_Kt[:,t,f])))
            if isnan(single_mean)==True: 
                p_tdiff_Kt_data[t,f] = 0.
                p_tdstd_Kt_data[t,f] = 0.
            else: 
                p_tdiff_Kt_data[t,f] = single_mean
                p_tdstd_Kt_data[t,f] = single_std


pylab.rc("font",size=9)
pylab.subplot(1,2,1)
pylab.imshow(p_tmean_Kt_data,vmin=0,vmax=8000,aspect=(new_freq[-1]-new_freq[0])/(times[-1]-times[0]),extent=(new_freq[0],new_freq[-1],times[-1],times[0]))
pylab.colorbar()
pylab.ylabel('LST (Hours)')
pylab.xlabel('Frequency (MHz)')
pylab.title('JNC Calibrated Avg Data')

pylab.subplot(1,2,2)
pylab.imshow(p_tmean_Kdgsm_data,vmin=0,vmax=8000,aspect=(new_freq[-1]-new_freq[0])/(times[-1]-times[0]),extent=(new_freq[0],new_freq[-1],times[-1],times[0]))
pylab.colorbar()
pylab.ylabel('LST (Hours)')
pylab.xlabel('Frequency (MHz)')
pylab.title('Delta GSM Calibrated Avg Data') 

pylab.savefig(outdir+'Combined_avg_mean_spectra',dpi=300)
pylab.clf()

pylab.rc("font",size=7)
for t in range(0,15):
    t10 = where(times<=t*1.7)[0][-1]
    pylab.subplot(4,4,t+1)
    pylab.scatter(new_freq,p_tmean_Kt_data[t10]/1000.,c='b',edgecolor='b',s=3)
    pylab.errorbar(new_freq,p_tmean_Kt_data[t10]/1000.,p_tstd_Kt_data[t10]/1000.,ecolor='b',fmt=None)
    pylab.scatter(new_freq,p_tmean_Kdgsm_data[t10]/1000.,c='g',edgecolor='g',s=3)
    pylab.errorbar(new_freq,p_tmean_Kdgsm_data[t10]/1000.,p_tstd_Kdgsm_data[t10]/1000.,ecolor='g',fmt=None)
    pylab.ylim(0,10)
    pylab.grid()
    pylab.xlim(50,100)
    if (t+1)%4 ==1:
        pylab.ylabel('Temperature (1000 K)')
    if t>11:
        pylab.xlabel('Frequency (MHz)')
    
    pylab.text(70,7,str(int(times[t10]))+'Hours LST',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))

pylab.subplot(4,4,16)
pylab.axis('off')
pylab.text(0,0.3,'Blue is Kt',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
pylab.text(0,0.7,'Green is Kdgsm',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
 
pylab.savefig(outdir+'Combined_K_frequency_spectra',dpi=300)
pylab.clf()

pylab.rc("font",size=7) 
for t in range(0,15): 
    t10 = where(times<=t*1.7)[0][-1]
    pylab.subplot(4,4,t+1) 
    pylab.scatter(new_freq,p_tdiff_Kt_data[t10],c='b',edgecolor='b',s=3)
    pylab.errorbar(new_freq,p_tdiff_Kt_data[t10],p_tdstd_Kt_data[t10],ecolor='b',fmt=None)
    pylab.scatter(new_freq,p_tdiff_Kdgsm_data[t10],c='g',edgecolor='g',s=3)
    pylab.errorbar(new_freq,p_tdiff_Kdgsm_data[t10],p_tdstd_Kdgsm_data[t10],ecolor='g',fmt=None)
    pylab.ylim(-100,100)
    pylab.grid()
    pylab.xlim(60,90)
    if (t+1)%4 ==1:
        pylab.ylabel('Residual T (K)')
    if t>11:
        pylab.xlabel('Frequency (MHz)')
 
    pylab.text(70,70,str(int(times[t10]))+'Hours LST',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
 
pylab.subplot(4,4,16)
pylab.axis('off')
pylab.text(0,0.3,'Blue is Kt',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
pylab.text(0,0.7,'Green is Kdgsm',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
 
pylab.savefig(outdir+'Combined_K_frequency_residuals',dpi=300)
pylab.clf()

#pylab.savefig(outdir+'Combined_K_frequency_spectra_test',dpi=300)
#pylab.clf()


pylab.rc("font",size=7)
freq_start=where(new_freq<=60.)[0][-1]
for i in range(freq_start,freq_start+15):
    pylab.subplot(4,4,i+1-freq_start)
    pylab.scatter(times,p_tmean_Kt_data[:,i]/1000.,c='b',edgecolor='b',s=1)
#    pylab.errorbar(times,p_tmean_Kt_data[:,i]/1000.,p_tstd_Kt_data[:,i]/1000.,c='b')
    pylab.scatter(times,p_tmean_Kdgsm_data[:,i]/1000.,c='g',edgecolor='g',s=1)
#    pylab.errorbar(times,p_tmean_Kdgsm_data[:,i]/1000.,p_tstd_Kdgsm_data[:,i]/1000.,c='g')
    pylab.xlim(0,24)
    
    if new_freq[i]<=57:
        pylab.ylim(0,10)
        pylab.text(1,8.8,str(int(new_freq[i]))+' MHz',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
    elif new_freq[i]<=65:
        pylab.ylim(0,8)
        pylab.text(1,7,str(int(new_freq[i]))+' MHz',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
    elif new_freq[i]<=75.:
        pylab.ylim(0,6)
        pylab.text(1,5.2,str(int(new_freq[i]))+' MHz',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
    else:
        pylab.ylim(0,4)
        pylab.text(1,3.4,str(int(new_freq[i]))+' MHz',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))

    pylab.grid()
    if (i+1-freq_start)%4 ==1:
        pylab.ylabel('Temperature (1000 K)')
    if i>freq_start+10:
        pylab.xlabel('LST (Hours)')

pylab.subplot(4,4,16)
pylab.axis('off')
pylab.text(0,0.3,'Blue is Kt',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
pylab.text(0,0.7,'Green is Kdgsm',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))

pylab.savefig(outdir+'Combined_K_time_spectra',dpi=300)
pylab.clf()


p_tot_mean_Kt = zeros(len(new_freq))
p_tot_std_Kt = zeros(len(new_freq))
p_tot_mean_Kdgsm = zeros(len(new_freq))
p_tot_std_Kdgsm = zeros(len(new_freq))
p_tot_diff_Kt = zeros(len(new_freq))
p_tot_diff_Kdgsm = zeros(len(new_freq))
p_tot_dstd_Kt = zeros(len(new_freq))
p_tot_dstd_Kdgsm = zeros(len(new_freq))
tmin = where(times<=13.)[0][-1]
tmax= where(times<=19.)[0][-1]
for f in range(0,len(new_freq)):
    p_tot_mean_Kt[f] = ma.mean(ma.compressed(ma.array(p_tmean_Kt_data[tmin:tmax,f],mask=p_tmask_Kt[tmin:tmax,f])))
    p_tot_mean_Kdgsm[f] = ma.mean(ma.compressed(ma.array(p_tmean_Kdgsm_data[tmin:tmax,f],mask=p_tmask_Kdgsm[tmin:tmax,f])))
    p_tot_std_Kt[f] = ma.std(ma.compressed(ma.array(p_tmean_Kt_data[tmin:tmax,f],mask=p_tmask_Kt[tmin:tmax,f])))
    p_tot_std_Kdgsm[f] = ma.std(ma.compressed(ma.array(p_tmean_Kdgsm_data[tmin:tmax,f],mask=p_tmask_Kdgsm[tmin:tmax,f])))
    p_tot_diff_Kt[f] = ma.mean(ma.compressed(ma.array(p_tdiff_Kt_data[tmin:tmax,f],mask=p_tmask_Kt[tmin:tmax,f])))
    p_tot_diff_Kdgsm[f] = ma.mean(ma.compressed(ma.array(p_tdiff_Kdgsm_data[tmin:tmax,f],mask=p_tmask_Kdgsm[tmin:tmax,f])))
    p_tot_dstd_Kt[f] = ma.std(ma.compressed(ma.array(p_tdiff_Kt_data[tmin:tmax,f],mask=p_tmask_Kt[tmin:tmax,f])))
    p_tot_dstd_Kdgsm[f] = ma.std(ma.compressed(ma.array(p_tdiff_Kdgsm_data[tmin:tmax,f],mask=p_tmask_Kdgsm[tmin:tmax,f])))

p_tot_fit_Kt,p_tot_params_Kt = cf.poly_fore(p_tot_mean_Kt,zeros(len(p_tot_mean_Kt)),new_freq,60.,90.,2,p_tot_std_Kt)
p_tot_fit_Kdgsm,p_tot_params_Kt = cf.poly_fore(p_tot_mean_Kdgsm,zeros(len(p_tot_mean_Kdgsm)),new_freq,60.,90.,2,p_tot_std_Kdgsm)

pylab.scatter(new_freq,p_tot_mean_Kt,c='b',edgecolor='b',s=5)
pylab.errorbar(new_freq,p_tot_mean_Kt,p_tot_std_Kt,ecolor='b',fmt=None)
pylab.scatter(new_freq,p_tot_mean_Kdgsm,c='g',edgecolor='g',s=5)
pylab.errorbar(new_freq,p_tot_mean_Kdgsm,p_tot_std_Kdgsm,ecolor='g',fmt=None)
pylab.plot(new_freq,p_tot_fit_Kt,'c')
pylab.plot(new_freq,p_tot_fit_Kdgsm,'r')
pylab.ylim(0,8000)
pylab.ylabel('Temperature (K)')
pylab.xlim(50,100)
pylab.xlabel('Frequency (MHz)')
pylab.grid()
pylab.text(70,7000,'Blue is Kt, Green is Kdgsm',fontsize=12,bbox=dict(edgecolor='black',facecolor='white'))
pylab.savefig(outdir+'Combined_avg_K_frequency_spectra',dpi=300)
pylab.clf()

pylab.scatter(new_freq,p_tot_diff_Kt,c='b',edgecolor='b',s=8)
pylab.errorbar(new_freq,p_tot_diff_Kt,p_tot_dstd_Kt,ecolor='b',fmt=None)
pylab.scatter(new_freq,p_tot_diff_Kdgsm,c='g',edgecolor='g',s=8)
pylab.errorbar(new_freq,p_tot_diff_Kdgsm,p_tot_dstd_Kdgsm,ecolor='g',fmt=None)
pylab.scatter(new_freq,p_tot_fit_Kt - p_tot_mean_Kt,c='c',edgecolor='c',s=8)
pylab.scatter(new_freq,p_tot_fit_Kdgsm-p_tot_mean_Kdgsm,c='r',edgecolor='r',s=8)
#pylab.errorbar(new_freq,p_tot_fit_Kt-p_tot_mean_Kt,p_tot_std_Kt,ecolor='c',fmt=None)
#pylab.errorbar(new_freq,p_tot_fit_Kdgsm-p_tot_mean_Kdgsm,p_tot_std_Kdgsm,ecolor='r',fmt=None)

pylab.ylim(-100,100)
pylab.ylabel('Residual T (K)')
pylab.xlim(60,90)
pylab.xlabel('Frequency (MHz)')
pylab.grid()
pylab.text(80,90,'Blue is Kt, Green is Kdgsm',fontsize=12,bbox=dict(edgecolor='black',facecolor='white'))
pylab.savefig(outdir+'Combined_avg_K_frequency_resid',dpi=300)
pylab.clf()

print 'Kt mean data: ',p_tot_mean_Kt
print 'Kdgsm mean data: ',p_tot_mean_Kdgsm
    
#fig = pylab.figure(figsize=(8,6),dpi=300)
#pylab.rc("font",size=6)

#processed_data = psingle_Kdgsm_data
f70 = where(freqs<=70.)[0][-1]
processed_data = pfull_Kt_data[:,:,f70]/1000.
print freqs[f70]
print shape(processed_data)
processed_time = []
for i in range(0,len(date_ind)):
    processed_time.append(times)

#for f in range(0,len(freqs)):
#processed_data = pfull_Kdgsm_data[:,:,t]/1000.
pylab.rc("font",size=7) 
for i in range(0,len(processed_data)):
    pylab.subplot(4,4,i+1)
    pylab.scatter(processed_time[i],processed_data[i],c='g',edgecolor='g',s=3)
    pylab.xlim(0,24) 
    pylab.ylim(0,8) 
    pylab.grid()
    if (i+1)%4==1:
        pylab.ylabel('Temperature (1000 K)') 
    if len(processed_data)-i <=4:
        pylab.xlabel('LST (Hours)')
    date_label='June '+str(int(i+1))
    pylab.legend([date_label],loc=2)

pylab.subplot(4,4,16)
freq_data = 'Frequency is '+str(int(freqs[f70]*10)/10.)+' MHz'
pylab.axis('off')
pylab.text(0,0.7,freq_data,fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
pylab.text(0,0.3,'Calibration is JNC',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
#pylab.savefig(outdir+'Combined_Kdgsm_time_series',dpi=300) 
pylab.savefig(outdir+'Combined_Kt_test_time_series',dpi=300)
pylab.clf() 

processed_data = pfull_Kdgsm_data[:,:,f70]/1000.
for i in range(0,len(processed_data)):
    pylab.subplot(4,4,i+1)
    pylab.scatter(processed_time[i],processed_data[i],c='g',edgecolor='g',s=3)
    pylab.xlim(0,24)
    pylab.ylim(0,8)
    pylab.grid()
    if (i+1)%4==1:
        pylab.ylabel('Temperature (1000 K)')
    if len(processed_data)-i <=4:
        pylab.xlabel('LST (Hours)')
    date_label='June '+str(int(i+1))
    pylab.legend([date_label],loc=2)

pylab.subplot(4,4,16)
freq_data = 'Frequency is '+str(int(freqs[f70]*10)/10.)+' MHz'
pylab.axis('off')
pylab.text(0,0.7,freq_data,fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
pylab.text(0,0.3,'Calibration is Delta GSM',fontsize=9,bbox=dict(edgecolor='black',facecolor='white'))
#pylab.savefig(outdir+'Combined_Kdgsm_time_series',dpi=300) 
pylab.savefig(outdir+'Combined_Kdgsm_test_time_series',dpi=300)
pylab.clf()

