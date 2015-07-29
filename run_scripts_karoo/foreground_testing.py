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

data = numpy.load(indir+'gsm_cal_data_masked_Apr_03_70MHz_ant.npy')
times = numpy.load(indir+'gsm_cal_times_Apr_03_70MHz_ant.npy')
Kdgsm = numpy.load(indir+'gsm_cal_values_Apr_03_70MHz_ant.npy')
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
st,sf,short_data,sm,sf,sv,ste = ff.loadsingle(supdir+'2015-04-05-00-06-26_Ch_2_noisetrunc.dat')

#Best fit parameters to short data collected in the Karoo.
fitshort = [9.43299312e-9,-1.16197388e-10,4.31005321e-13]
short_fit = poly.polyval(freqs,fitshort)

smask = zeros(len(freqs))
rebin_short,short_mask,short_f = ff.rebin((short_data)*Kdgsm,smask,freqs,fscale)
smooth_short = itp.UnivariateSpline(freqs,short_data*Kdgsm)
short_sm = smooth_short(freqs)
new_data = data - short_sm 
srdata = rebin_data-smooth_short(short_f)
f50 = where(freqs<=50.)[0][-1]
f90 = where(freqs<=90.)[0][-1]
rf50 = where(short_f<=50.)[0][-1]
rf90 = where(short_f<=90.)[0][-1]
ns = [2,3] 
ms = [1,2,3,4] 
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

for r in ns:
    Kparams= zeros((len(data),r))
    sKparams=zeros((len(data),r))
    Rparams= zeros((len(data),len(ms),r+ms[-1]+1))
    sRparams=zeros((len(data),len(ms),r+ms[-1]+1))

    for t in range(0,len(times)):
        if len(mask[t])!=sum(mask[t]):
            Kfit[t], Kparams[t] = cf.poly_fore(data[t],mask[t],freqs,50.,90.,r-1,ones(len(data[t])))
            sKfit[t], sKparams[t] = cf.poly_fore(new_data[t],mask[t],freqs,50.,90.,r-1,ones(len(data[t])))
            for m in ms:
                Rfit[t,:,r-ns[0],m-ms[0]], Rp = cf.rat_fore(data[t],mask[t],freqs,50.,90.,r-1,m)
#                print shape(Rp)
                Rparams[t,m-ms[0],0:(len(Rp))] = Rp
                sRfit[t,:,r-ns[0],m-ms[0]], sRp = cf.rat_fore(new_data[t],mask[t],freqs,50.,90.,r-1,m)
#                print shape(sRp)
                sRparams[t,m-ms[0],0:(len(sRp))] = sRp
        rdata[t],rmask[t],nf = ff.rebin(data[t]-Kfit[t],mask[t],freqs,fscale)
        nandata = where(isnan(rdata))
        rdata[nandata]=0.
        srdata[t],srmask[t],nf = ff.rebin(new_data[t]-sKfit[t],mask[t],freqs,fscale)
        nandata = where(isnan(srdata))
        srdata[nandata] = 0.
        for m in ms:
            rrdata[t,:,m-ms[0]],rrmask[t,:,m-ms[0]],nf = ff.rebin(data[t]-Rfit[t,:,r-ns[0],m-ms[0]],mask[t],freqs,fscale) 
            rsrdata[t,:,m-ms[0]],rsrmask[t,:,m-ms[0]],nf = ff.rebin(new_data[t]-sRfit[t,:,r-ns[0],m-ms[0]],mask[t],freqs,fscale)
    
    nandata = where(isnan(rrdata))
    rrdata[nandata] = 0.0
    nandata = where(isnan(rsrdata))
    rsrdata[nandata] = 0.0

    mean_Kf[:,r-ns[0]],mean_Kmask[:,r-ns[0]] = cf.time_mean(Kfit,mask)
    smean_Kf[:,r-ns[0]],smean_Kmask[:,r-ns[0]] = cf.time_mean(sKfit,mask)
    mean_Kresid[r-ns[0],:],mean_Krmask[r-ns[0],:] = cf.time_mean(rdata,rmask)
    mean_sKresid[r-ns[0],:],mean_sKrmask[r-ns[0],:] = cf.time_mean(srdata,srmask)

    for i in range(0,len(ms)):
        mean_Rf[:,r-ns[0],i],mean_Rmask[:,r-ns[0],i] = cf.time_mean(Rfit[:,:,r-ns[0],i],mask)
        mean_Rresid[r-ns[0],i,:],mean_Rrmask[r-ns[0],i,:] = cf.time_mean(rrdata[:,:,i],rrmask) 
        smean_Rf[:,r-ns[0],i],smean_Rmask[:,r-ns[0],i] = cf.time_mean(sRfit[:,:,r-ns[0],i],mask)
        mean_sRresid[r-ns[0],i,:],mean_sRrmask[r-ns[0],i,:] = cf.time_mean(rsrdata[:,:,i],rsrmask)

pylab.rc('font',size=8)
index = 1
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        pylab.subplot(len(ns),len(ms),index)
        pylab.scatter(ma.compressed(ma.array(freqs,mask=mean_smask)),ma.compressed(ma.array(mean_sig,mask=mean_smask)),s=1,c='b',edgecolor='b')
        pylab.plot(freqs,mean_Kf[:,n],c='c')
        pylab.plot(freqs,mean_Rf[:,n,m],c='g')

        pylab.scatter(ma.compressed(ma.array(freqs,mask=smean_smask)),ma.compressed(ma.array(smean_sig,mask=smean_smask)),s=1,c='r',edgecolor='r') 
        pylab.plot(freqs,smean_Kf[:,n],c='k')
        pylab.plot(freqs,smean_Rf[:,n,m],c='m') 
#        pylab.scatter(ma.compressed(ma.array(freqs,mask=mean_mask)),ma.compressed(ma.array(short_sm,mask=mean_mask)),s=1,c='y',edgecolor='y')
        pylab.xlim(60,90)
        pylab.ylim(1e3,4e3)
        pylab.grid()
        if (index)%len(ms)==1:
            pylab.ylabel('Temperature (Kelvin)')
        if (index)>((len(ns)-1)*len(ms)):
            pylab.xlabel('Frequency (MHz)')
        index+=1
pylab.savefig(outdir+'total_fit_data.png',dpi=300)
pylab.clf()

index=1
pylab.rc('font',size=8)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        pylab.subplot(len(ns),len(ms),index)
        pylab.scatter(ma.compressed(ma.array(nf,mask=mean_Krmask[n,:])),ma.compressed(ma.array(mean_Kresid[n,:],mask=mean_Krmask[n,:])),s=1,c='b',edgecolor='b')
        pylab.scatter(ma.compressed(ma.array(nf,mask=mean_sKrmask[n,:])),ma.compressed(ma.array(mean_sKresid[n,:],mask=mean_sKrmask[n,:])),s=1,c='c',edgecolor='c')
        pylab.scatter(ma.compressed(ma.array(nf,mask=mean_Rrmask[n,m])),ma.compressed(ma.array(mean_Rresid[n,m],mask=mean_Rrmask[n,m])),s=1,c='r',edgecolor='r')
        pylab.scatter(ma.compressed(ma.array(nf,mask=mean_sRrmask[n,m])),ma.compressed(ma.array(mean_sRresid[n,m],mask=mean_sRrmask[n,m])),s=1,c='m',edgecolor='m')
        mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
        std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
        pylab.xlim(50,90)
        pylab.ylim(mean_val-std_val,mean_val+std_val)
        pylab.grid()

        if (index)%len(ms)==1:
            pylab.ylabel('Resid Temp (Kelvin)')
        if (index)>(len(ns)-1)*len(ms):
            pylab.xlabel('Frequency (MHz)')
        index+=1
pylab.savefig(outdir+'total_resid_data.png',dpi=300)
pylab.clf()

index=1
pylab.rc('font',size=8)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        pylab.subplot(len(ns),len(ms),index)
        mean_val = ma.mean(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
        std_val = ma.std(ma.compressed(ma.array(mean_Kresid[n,rf50:rf90],mask=mean_Krmask[n,rf50:rf90])))
        pylab.imshow(data[:,f50:f90]-Rfit[:,f50:f90,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=40./24.,extent=(freqs[f50],freqs[f90],times[-1],times[0]))
        cb = pylab.colorbar()
        if (index)%len(ms)==(len(ms)):
            cb.set_label('Temp Resid (K)')
        if (index)%len(ms)==1:
            pylab.ylabel('Sidereal Time (Hrs)')
        if (index)>((len(ns)-1)*len(ms)):
            pylab.xlabel('Frequency (MHz)')
        index+=1
pylab.savefig(outdir+'R_resid_distribution_data.png',dpi=300) 
pylab.clf()

index=1
pylab.rc('font',size=8)
for n in range(0,len(ns)):
    for m in range(0,len(ms)):
        pylab.subplot(len(ns),len(ms),index)
        mean_val = ma.mean(ma.compressed(ma.array(mean_sKresid[n,rf50:rf90],mask=mean_sKrmask[n,rf50:rf90])))
        std_val = ma.std(ma.compressed(ma.array(mean_sKresid[n,rf50:rf90],mask=mean_sKrmask[n,rf50:rf90])))
        pylab.imshow(new_data[:,f50:f90]-sRfit[:,f50:f90,n,m],vmin=mean_val-std_val,vmax=mean_val+std_val,aspect=40./24.,extent=(freqs[f50],freqs[f90],times[-1],times[0]))
        cb = pylab.colorbar()
        if (index)%len(ns)==(len(ns)): 
            cb.set_label('Temp Resid (K)')
        if (index)%len(ns)==1: 
            pylab.ylabel('Sidereal Time (Hrs)')
        if (index)>((len(ns)-1)*len(ms)):
            pylab.xlabel('Frequency (MHz)')
        index+=1
pylab.savefig(outdir+'short_R_resid_distribution_data.png',dpi=300)
pylab.clf()


