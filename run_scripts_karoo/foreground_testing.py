"""
Module for testing out new code to improve foreground subtraction.
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
direct = sys.argv[3]
ant = sys.argv[3].split('_')[-1]
supdir= '../../supplemental_data/'

lat = '-30.727206'
lon = '21.479055'
elevation = 1080
idate = '2015/04/01'
site = eph.Observer()
site.lon = lon
site.lat = lat
site.elevation = elevation
date = eph.date(idate)
site.date = date
site.pressure = 0
sun = eph.Sun(site)
site.date = site.next_rising(sun)
print site.sidereal_time()
site.date = site.next_setting(sun)
print site.sidereal_time()

data = numpy.load(indir+'gsm_cal_data_masked_Apr_'+direct+'MHz_ant.npy')
times = numpy.load(indir+'gsm_cal_times_Apr_'+direct+'MHz_ant.npy')
Kdgsm = numpy.load(indir+'gsm_cal_values_Apr_'+direct+'MHz_ant.npy')
freqs = arange(40.,130.,90./len(data[0]))

mask = zeros((len(data),len(data[0])))
test = where(data==0.0)
mask[test] = 1.0
tmask=zeros(len(data))
for t in range(0,len(times)):
    if len(mask[t])==sum(mask[t]):
        tmask[t] = 1.
fscale = 16

rebin_data = zeros((len(data),len(data[0])/fscale))
rebin_mask = zeros((len(data),len(data[0])/fscale))
for t in range(0,len(times)):
    rebin_data[t],rebin_mask[t],rebin_f = ff.rebin(data[t],mask[t],freqs,fscale)
st,sf,short_data,sm,sfreq,sv,ste = ff.loadsingle(supdir+'2015-04-05-00-06-26_Ch_2_noisetrunc.dat')

#Best fit parameters to short data collected in the Karoo.
#fitshort = [9.43299312e-9,-1.16197388e-10,4.31005321e-13]
#short_fit = poly.polyval(freqs,fitshort)

fitshort = [4.47140321e-11,3.09422597,4.15995054,9.43146423e-09,-1.16198506e-10,4.31190367e-13]
def funci(x,a,b,c,d,e,f):
    return a*cos(2*pi*x/b+c)+d+e*x+f*x**2
fitshort,cov=opt.curve_fit(funci,freqs,short_data,fitshort[:])
short_fit = funci(freqs,fitshort[0],fitshort[1],fitshort[2],fitshort[3],fitshort[4],fitshort[5])

smask = zeros(len(freqs))
rebin_short,short_mask,short_f = ff.rebin((short_data)*Kdgsm,smask,freqs,fscale)
#smooth_short = itp.UnivariateSpline(freqs,short_data*Kdgsm)
#short_sm = smooth_short(freqs)
#new_data = data - short_sm 
new_data = data - short_fit*Kdgsm
srdata = rebin_data-funci(short_f,fitshort[0],fitshort[1],fitshort[2],fitshort[3],fitshort[4],fitshort[5])
f50 = where(freqs<=50.)[0][-1]
f90 = where(freqs<=90.)[0][-1]
f80 = where(freqs<=80.)[0][-1]
f110 = where(freqs<=110.)[0][-1]
rf50 = where(short_f<=50.)[0][-1]
rf90 = where(short_f<=90.)[0][-1]
rf80 = where(short_f<=80.)[0][-1]
rf110 = where(short_f<=110.)[0][-1]
ns = [2,3,4] 
ms = [1,2,3] 
fscale = 16

mean_sig,mean_smask = cf.time_mean(data,mask)
rmean_sig,rmean_mask = cf.time_mean(rebin_data,rebin_mask)
smean_sig,smean_smask = cf.time_mean(new_data,mask)
srmean_sig,srmean_mask = cf.time_mean(srdata,rebin_mask)

Kfit = zeros((len(data),len(data[0])))
sKfit = zeros((len(data),len(data[0])))
Rfit = zeros((len(data),len(data[0]),len(ns),len(ms)))
sRfit = zeros((len(data),len(data[0]),len(ns),len(ms)))

rdata = zeros((len(data),len(freqs)/fscale))
rmask = zeros((len(data),len(freqs)/fscale))
srdata = zeros((len(data),len(freqs)/fscale))
srmask = zeros((len(data),len(freqs)/fscale))
rrdata = zeros((len(data),len(freqs)/fscale,len(ms)))
rrmask = zeros((len(data),len(freqs)/fscale,len(ms)))
rsrdata = zeros((len(data),len(freqs)/fscale,len(ms)))
rsrmask = zeros((len(data),len(freqs)/fscale,len(ms)))

mean_Kresid = zeros((len(ns),len(freqs)/fscale))
mean_Krmask = zeros((len(ns),len(freqs)/fscale))
mean_sKresid = zeros((len(ns),len(freqs)/fscale))
mean_sKrmask = zeros((len(ns),len(freqs)/fscale))
mean_Rresid = zeros((len(ns),len(ms),len(freqs)/fscale))
mean_Rrmask = zeros((len(ns),len(ms),len(freqs)/fscale))
mean_sRresid = zeros((len(ns),len(ms),len(freqs)/fscale))
mean_sRrmask = zeros((len(ns),len(ms),len(freqs)/fscale))

mean_Kf = zeros((len(freqs),len(ns)))
mean_Kmask = zeros((len(freqs),len(ns)))
smean_Kf = zeros((len(freqs),len(ns)))
smean_Kmask = zeros((len(freqs),len(ns)))
mean_Rf = zeros((len(freqs),len(ns),len(ms)))
mean_Rmask = zeros((len(freqs),len(ns),len(ms)))
smean_Rf = zeros((len(freqs),len(ns),len(ms)))
smean_Rmask = zeros((len(freqs),len(ns),len(ms)))

fRparams = zeros((len(data),len(ns),len(ms),5))
fRdata = zeros((len(data),len(freqs),len(ns),len(ms)))

def func(x,a,b,c,d,e):
    return abs(a)*cos(2*pi*x/b+c%(2*pi))+d+e*x
rf_params = [5.,3.,0.01,10.,1.]
fit_min = 50
fit_max = 110
for r in ns:
    Kparams= zeros((len(data),r))
    sKparams=zeros((len(data),r))
    Rparams= zeros((len(data),len(ms),r+ms[-1]+1))
    sRparams=zeros((len(data),len(ms),r+ms[-1]+1))

    for t in range(0,len(times)):
        if len(mask[t])!=sum(mask[t]):
            if ant=='70':
                fit_min = 50.
                fit_max = 90.
            elif ant=='100':
                fit_min = 80.
                fit_max = 110.

            Kfit[t], Kparams[t] = cf.poly_fore(data[t],mask[t],freqs,fit_min,fit_max,r-1,ones(len(data[t])))
            sKfit[t], sKparams[t] = cf.poly_fore(new_data[t],mask[t],freqs,fit_min,fit_max,r-1,ones(len(data[t])))
            for m in ms:
                Rfit[t,:,r-ns[0],m-ms[0]], Rp = cf.rat_fore(data[t],mask[t],freqs,fit_min,fit_max,r-1,m)
                Rparams[t,m-ms[0],0:(len(Rp))] = Rp
                sRfit[t,:,r-ns[0],m-ms[0]], sRp = cf.rat_fore(new_data[t],mask[t],freqs,fit_min,fit_max,r-1,m)
                sRparams[t,m-ms[0],0:(len(sRp))] = sRp
        rdata[t],rmask[t],nf = ff.rebin(data[t]-Kfit[t],mask[t],freqs,fscale)
        srdata[t],srmask[t],nf = ff.rebin(new_data[t]-sKfit[t],mask[t],freqs,fscale)
        rf_params = [5.,3.,0.01,10.,1.]

        for m in ms:
            rrdata[t,:,m-ms[0]],rrmask[t,:,m-ms[0]],nf = ff.rebin(data[t]-Rfit[t,:,r-ns[0],m-ms[0]],mask[t],freqs,fscale) 
            rsrdata[t,:,m-ms[0]],rsrmask[t,:,m-ms[0]],nf = ff.rebin(new_data[t]-sRfit[t,:,r-ns[0],m-ms[0]],mask[t],freqs,fscale)
 
            rfmin = where(freqs<=(fit_min+15.))[0][-1]
            rfmax = where(freqs<=(fit_min+20.))[0][-1]
            if ma.mean(data[t,rfmin:rfmax])>0.0:
                try:
                    rf_params,cov = opt.curve_fit(func,freqs[rfmin:rfmax],data[t,rfmin:rfmax]-Rfit[t,rfmin:rfmax,r-ns[0],m-ms[0]],rf_params[:],maxfev=5000)

                    fRparams[t,r-ns[0],m-ms[0]] = rf_params
                    fRdata[t,:,r-ns[0],m-ms[0]] = abs(rf_params[0])*cos(2*pi*freqs/rf_params[1]+rf_params[2]%(2*pi))
                except RuntimeError:
                    print 'fit failed at time ',t,' m value ',m,' and n value ',r
                
#        print t,r,m
            
    print ma.mean(fRparams,axis=0)
    print ma.std(fRparams,axis=0)
    mean_Kf[:,r-ns[0]],mean_Kmask[:,r-ns[0]] = cf.time_mean(Kfit,mask)
    smean_Kf[:,r-ns[0]],smean_Kmask[:,r-ns[0]] = cf.time_mean(sKfit,mask)
    mean_Kresid[r-ns[0],:],mean_Krmask[r-ns[0],:] = cf.time_mean(rdata,rmask)
    mean_sKresid[r-ns[0],:],mean_sKrmask[r-ns[0],:] = cf.time_mean(srdata,srmask)

    for i in range(0,len(ms)):
        mean_Rf[:,r-ns[0],i],mean_Rmask[:,r-ns[0],i] = cf.time_mean(Rfit[:,:,r-ns[0],i],mask)
        mean_Rresid[r-ns[0],i,:],mean_Rrmask[r-ns[0],i,:] = cf.time_mean(rrdata[:,:,i],rrmask[:,:,i]) 
        smean_Rf[:,r-ns[0],i],smean_Rmask[:,r-ns[0],i] = cf.time_mean(sRfit[:,:,r-ns[0],i],mask)
        mean_sRresid[r-ns[0],i,:],mean_sRrmask[r-ns[0],i,:] = cf.time_mean(rsrdata[:,:,i],rsrmask[:,:,i])


index=1
pylab.rc('font',size=6)
f,axarr = plt.subplots(len(ns),len(ms),sharex='col',sharey='row')
f.suptitle('Time Mean Data and Polynomial Fits, n=num order, m=denom order',size=10)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        axarr[n,m].scatter(ma.compressed(ma.array(freqs,mask=mean_smask)),ma.compressed(ma.array(mean_sig,mask=mean_smask))/1000.,s=1,c='b',edgecolor='b',label='data')
        axarr[n,m].scatter(ma.compressed(ma.array(freqs,mask=mean_smask)),ma.compressed(ma.array(mean_Kf[:,n]/1000.,mask=mean_smask)),c='c',edgecolor='c',s=1,label='P')
        axarr[n,m].scatter(ma.compressed(ma.array(freqs,mask=mean_smask)),ma.compressed(ma.array(mean_Rf[:,n,m]/1000.,mask=mean_smask)),c='g',edgecolor='g',s=1,label='R')

        axarr[n,m].scatter(ma.compressed(ma.array(freqs,mask=smean_smask)),ma.compressed(ma.array(smean_sig,mask=smean_smask))/1000.,s=1,c='r',edgecolor='r',label='sh data') 
        axarr[n,m].scatter(ma.compressed(ma.array(freqs,mask=smean_smask)),ma.compressed(ma.array(smean_Kf[:,n]/1000.,mask=smean_smask)),c='k', edgecolor='k',s=1,label='sh P')
        axarr[n,m].scatter(ma.compressed(ma.array(freqs,mask=smean_smask)),ma.compressed(ma.array(smean_Rf[:,n,m]/1000.,mask=smean_smask)),c='m',edgecolor='m',s=1,label='sh R') 
#        pylab.scatter(ma.compressed(ma.array(freqs,mask=mean_mask)),ma.compressed(ma.array(short_sm,mask=mean_mask)),s=1,c='y',edgecolor='y')
        axarr[n,m].grid()
        axarr[n,m].legend(prop={'size':6})
        axarr[n,m].set_title('n='+str(ns[n]-1)+', m='+str(ms[m]))
        if (index)%len(ms)==1:
            axarr[n,m].set_ylabel('Temp (1000 Kelvin)')
            axarr[n,m].set_ylim(0,6)
        if (index)>((len(ns)-1)*len(ms)):
            axarr[n,m].set_xlabel('Frequency (MHz)')
            if ant=='100':
                axarr[n,m].set_xlim(70,120)
            elif ant=='70':
                axarr[n,m].set_xlim(40,100)
        index+=1
f.subplots_adjust(wspace=0)
pylab.savefig(outdir+'total_fit_data.png',dpi=300)
pylab.clf()

index=1
pylab.rc('font',size=6)
f,axarr = plt.subplots(len(ns),len(ms),sharex='col',sharey='row')
f.suptitle('Time Mean Residuals, n=num order,m=denom order',size=10)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        axarr[n,m].scatter(ma.compressed(ma.array(nf,mask=mean_Krmask[n,:])),ma.compressed(ma.array(mean_Kresid[n,:],mask=mean_Krmask[n,:])),s=1,c='b',edgecolor='b',label='P')
        axarr[n,m].scatter(ma.compressed(ma.array(nf,mask=mean_sKrmask[n,:])),ma.compressed(ma.array(mean_sKresid[n,:],mask=mean_sKrmask[n,:])),s=1,c='c',edgecolor='c',label='sh P')
        axarr[n,m].scatter(ma.compressed(ma.array(nf,mask=mean_Rrmask[n,m])),ma.compressed(ma.array(mean_Rresid[n,m],mask=mean_Rrmask[n,m])),s=1,c='r',edgecolor='r',label='R')
        axarr[n,m].scatter(ma.compressed(ma.array(nf,mask=mean_sRrmask[n,m])),ma.compressed(ma.array(mean_sRresid[n,m],mask=mean_sRrmask[n,m])),s=1,c='m',edgecolor='m',label='sh R')
        if ant=='100':
            mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf80:rf110],mask=mean_Krmask[n,rf80:rf110])))
            std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf80:rf110],mask=mean_Krmask[n,rf80:rf110])))
        elif ant=='70':
            mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
            std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))

        axarr[n,m].grid() 
        axarr[n,m].legend(prop={'size':6})
        axarr[n,m].set_title('n='+str(ns[n]-1)+', m='+str(ms[m]))
        if (index)%len(ms)==1:
            axarr[n,m].set_ylabel('Resid Temp (Kelvin)')
            axarr[n,m].set_ylim(mean_val-std_val,mean_val+std_val)
        if (index)>(len(ns)-1)*len(ms):
            axarr[n,m].set_xlabel('Frequency (MHz)')
            if ant=='100':
                axarr[n,m].set_xlim(70,120)
            elif ant=='70': 
                axarr[n,m].set_xlim(40,100)
        index+=1
f.subplots_adjust(wspace=0)
pylab.savefig(outdir+'total_resid_data.png',dpi=300)
pylab.clf()

index=1
pylab.rc('font',size=7)
f,axarr = plt.subplots(len(ns),len(ms),sharex='col',sharey='row')
f.suptitle('Rational Function Residuals, n=num order,m=denom order',size=10)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        if ant=='100':
            mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf80:rf110],mask=mean_Krmask[n,rf80:rf110])))
            std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf80:rf110],mask=mean_Krmask[n,rf80:rf110])))
            im = axarr[n,m].imshow(data[:,f80:f110]-Rfit[:,f80:f110,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=30./24.,extent=(freqs[f80],freqs[f110],times[-1],times[0]))
        elif ant=='70':
            mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
            std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
            im = axarr[n,m].imshow(data[:,f50:f90]-Rfit[:,f50:f90,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=40./24.,extent=(freqs[f50],freqs[f90],times[-1],times[0]))
        axarr[n,m].set_title('n='+str(ns[n]-1)+', m='+str(ms[m])) 
        if (index)%len(ms)==0:
            cb = plt.colorbar(im,ax=axarr[n,m])
            cb.set_label('Temp Resid (K)')
        if (index)%len(ms)==1:
            axarr[n,m].set_ylabel('Sidereal Time (Hrs)')
        if (index)>((len(ns)-1)*len(ms)):
            axarr[n,m].set_xlabel('Frequency (MHz)')
        index+=1
f.subplots_adjust(wspace=0)
pylab.savefig(outdir+'R_resid_distribution_data.png',dpi=300) 
pylab.clf()

index=1
pylab.rc('font',size=7)
f,axarr = plt.subplots(len(ns),len(ms),sharex='col',sharey='row')
f.suptitle('Short Sub Rational Function Residuals, n=num order, m=denom order',size=10)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        if ant=='100':
            mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf80:rf110],mask=mean_Krmask[n,rf80:rf110])))
            std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf80:rf110],mask=mean_Krmask[n,rf80:rf110])))
            im = axarr[n,m].imshow(new_data[:,f80:f110]-sRfit[:,f80:f110,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=30./24.,extent=(freqs[f80],freqs[f110],times[-1],times[0]))

        elif ant=='70':
            mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
            std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
            im = axarr[n,m].imshow(new_data[:,f50:f90]-sRfit[:,f50:f90,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=40./24.,extent=(freqs[f50],freqs[f90],times[-1],times[0]))
        axarr[n,m].set_title('n='+str(ns[n]-1)+', m='+str(ms[m])) 
        if (index)%len(ms)==0: 
            cb = plt.colorbar(im,ax=axarr[n,m])
            cb.set_label('Temp Resid (K)')
        if (index)%len(ms)==1: 
            axarr[n,m].set_ylabel('Sidereal Time (Hrs)')
        if (index)>((len(ns)-1)*len(ms)):
            axarr[n,m].set_xlabel('Frequency (MHz)')
        index+=1
f.subplots_adjust(wspace=0) 
pylab.savefig(outdir+'short_R_resid_distribution_data.png',dpi=300)
pylab.clf()

index=1
pylab.rc('font',size=7) 
f,axarr = plt.subplots(len(ns),len(ms),sharex='col',sharey='row')
f.suptitle('Rational Function Resid w/ Sin remove, n=num order, m=denom order',size=10)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        if ant=='100':
            mean_val = ma.mean(ma.compressed(ma.array(data[:,f80:f110]-Rfit[:,f80:f110,n,m],mask=mask[:,f80:f110])))
            std_val = ma.std(ma.compressed(ma.array(data[:,f80:f110]-Rfit[:,f80:f110,n,m],mask=mask[:,f80:f110])))
            im = axarr[n,m].imshow((data[:,f80:f110]-Rfit[:,f80:f110,n,m])-fRdata[:,f80:f110,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=30./24.,extent=(freqs[f80],freqs[f110],times[-1],times[0]))

        elif ant=='70': 
            mean_val = ma.mean(ma.compressed(ma.array(data[:,f50:f90]-Rfit[:,f50:f90,n,m],mask=mask[:,f50:f90])))
            std_val = ma.std(ma.compressed(ma.array(data[:,f50:f90]-Rfit[:,f50:f90,n,m],mask=mask[:,f50:f90])))
            im = axarr[n,m].imshow((data[:,f50:f90]-Rfit[:,f50:f90,n,m])-fRdata[:,f50:f90,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=40./24.,extent=(freqs[f50],freqs[f90],times[-1],times[0]))
        axarr[n,m].set_title('n='+str(ns[n]-1)+', m='+str(ms[m]))
        if (index)%len(ms)==0:
            cb = plt.colorbar(im,ax=axarr[n,m]) 
            cb.set_label('Temp Resid (K)')
        if (index)%len(ms)==1: 
            axarr[n,m].set_ylabel('Sidereal Time (Hrs)')
        if (index)>((len(ns)-1)*len(ms)):
            axarr[n,m].set_xlabel('Frequency (MHz)')
        index+=1
f.subplots_adjust(wspace=0) 
pylab.savefig(outdir+'R_resid_cos_fit_distribution_data.png',dpi=300)
pylab.clf()

index=1
pylab.rc('font',size=6)
f,axarr = plt.subplots(len(ns),len(ms),sharex='col',sharey='row')
f.suptitle('Cosine Fit Parameters, n=num order,m=denom order',size=10)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        axarr[n,m].scatter(ma.compressed(ma.array(times,mask=tmask)),ma.compressed(ma.array(fRparams[:,n,m,0],mask=tmask)),s=1,c='b',edgecolor='b',label='a (K) ')
        axarr[n,m].scatter(ma.compressed(ma.array(times,mask=tmask)),ma.compressed(ma.array(fRparams[:,n,m,1],mask=tmask)),s=1,c='c',edgecolor='c',label='b (MHz^-1)')
        axarr[n,m].scatter(ma.compressed(ma.array(times,mask=tmask)),ma.compressed(ma.array(fRparams[:,n,m,2]%(2*pi),mask=tmask)),s=1,c='r',edgecolor='r',label='c (radians)')
        axarr[n,m].grid()
        axarr[n,m].legend(prop={'size':6})
        axarr[n,m].set_title('n='+str(ns[n]-1)+', m='+str(ms[m]))
        if (index)%len(ms)==1:
            axarr[n,m].set_ylabel('Param Value')
            axarr[n,m].set_ylim(-20,20)
        if (index)>(len(ns)-1)*len(ms):
            axarr[n,m].set_xlabel('Sidereal Time (hrs)')
        index+=1
f.subplots_adjust(wspace=0)
pylab.savefig(outdir+'cos_param_data.png',dpi=300)
pylab.clf()

