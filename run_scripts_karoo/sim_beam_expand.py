"""
Code for making the initial loading files for gsm_generate_opt.py
Converts .dat files data into .npy arrays
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
import os
import sys
sys.path.append(os.path.abspath('../../hibiscus'))
import time

def load_sim_files(antenna,supdir):
    if antenna=='70':
        sim_data = np.load(supdir+'beam_simulations50-90.npy')
    elif antenna=='100':
        sim_data = np.load(supdir+'beam_simulations80-110.npy')

    gsm_freq = np.arange(50,111,1)

    curr_az = np.load(supdir+'sim_az.npy')
    curr_alt = np.load(supdir+'sim_alt.npy')

    return sim_data,curr_az,curr_alt,gsm_freq

def load_expand_files(antenna,indir):
    supdir = indir
    if antenna=='70':
        phi0 = np.loadtxt(supdir+'beam_sim_update/70MHz_Phi_0.dat')
        phi90 = np.loadtxt(supdir+'beam_sim_update/70MHz_Phi_90.dat')
    elif antenna=='100': 
        phi0 = np.loadtxt(supdir+'beam_sim_update/100MHz_Phi0.dat')
        phi90 = np.loadtxt(supdir+'beam_sim_update/100MHz_Phi90.dat')
 
    exp_alt = np.zeros(len(phi0)*2)
    exp_az = np.zeros(len(phi0)*2)
    exp_data = np.zeros((len(phi0)*2,len(phi0[0])-1))

    for i in range(0,len(phi0)): 
        if phi0[i,0]<=0:
            exp_alt[i] = np.pi/2.+phi0[i,0]*np.pi/180.
            exp_az[i] = np.pi
        else: 
            exp_alt[i] = np.pi/2.-phi0[i,0]*np.pi/180.
            exp_az[i] = 0.
        if phi90[i,0]<=0:
            exp_alt[i+len(phi0)] = np.pi/2.+phi0[i,0]*np.pi/180.
            exp_az[i+len(phi0)] = np.pi/2.
        else:
            exp_alt[i+len(phi0)] = np.pi/2.-phi90[i,0]*np.pi/180.
            exp_az[i+len(phi0)] = 3*np.pi/2
        exp_data[i] = 10**(phi0[i,1:]/10.)
        exp_data[i+len(phi0)] = 10**(phi90[i,1:]/10.)

    return exp_data,exp_az,exp_alt

def calc_sins(exp_data):
    offset = np.zeros((10,91))
    amp = np.zeros((10,91))
    for alt in range(0,91):
        for f in range(0,10):
            amp[f,alt] = abs(-exp_data[alt,f]+exp_data[alt+len(exp_data)/2,f])/2.
            if exp_data[alt,f]>=exp_data[alt+len(exp_data)/2,f]:
                offset[f,alt] = amp[f,alt]+exp_data[alt+len(exp_data)/2,f]
            else:
                offset[f,alt] = exp_data[alt,f]+amp[f,alt]
    return offset,amp

def calc_2nd_sin(antenna,sim_data,sim_alt,offset,amp):
    amp2 = np.zeros((len(amp),len(amp[0])))
    if antenna=='70':
        for i in range(0,91):
            inds = np.where(abs(sim_alt-i*np.pi/180.)<0.01)
            damp = (ma.max(10**(sim_data[40,inds]/10.))-ma.min(10**(sim_data[40,inds]/10.)))/2.
            doff = damp+ma.min(10**(sim_data[40,inds]/10.))    
            offset[7,i] = doff 
            amp[7,i] = damp
    elif antenna=='100':
        for i in range(0,91):
            inds = np.where(abs(sim_alt-i*np.pi/180.)<0.01)
            damp = (ma.max(10**(sim_data[30,inds]/10.))-ma.min(10**(sim_data[30,inds]/10.)))/2.
            doff = damp+ma.min(10**(sim_data[30,inds]/10.))    
            offset[6,i] = doff 
            amp[6,i] = damp

    fitfunc = lambda p,x: p[0]*np.exp(-(p[1]-x)**2/(2*p[2]**2))
    errfunc = lambda p,x,y: fitfunc(p,x)-y
    for i in range(0,91):
        for j in range(0,10):
            p,err = opt.leastsq(errfunc,[ma.max(amp[j]),np.where(amp[j]==ma.max(amp[j]))[0],10.],args=(np.arange(0,91),amp[j]))
            b = np.where(amp[j]==ma.max(amp[j]))[0]
            c = p[2]-3.
            amp2[j,i] = fitfunc([-(ma.max(amp[j])-1.1),b-9,c],i)/2.

    return offset, amp, amp2
 

def lin_freq_fit(antenna,offset,amp,amp2,sim_freq,sim_data,az):
    data = np.zeros(len(sim_freq))
    if antenna=='70':
        for f in range(0,len(sim_freq)):
            if sim_freq[f]<90.:
                data[f] = 10**(sim_data[f]/10.)
            elif sim_freq[f]==90.:
                of = offset[7]
                am = amp[7]
                am2 = amp2[7]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]<100.:
                m1 = (offset[8]-offset[7])/(sim_freq[50]-sim_freq[40])
                b1 = offset[7]-m1*sim_freq[40]
                of = m1*sim_freq[f]+b1
                m2 = (amp[8]-amp[7])/(sim_freq[50]-sim_freq[40])
                b2 = amp[7]-m2*sim_freq[40]
                am = m2*sim_freq[f]+b2
                m3 = (amp2[8]-amp2[7])/(sim_freq[50]-sim_freq[40]) 
                b3 = amp2[7]-m3*sim_freq[40]
                am2 = m3*sim_freq[f]+b3               
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]==100.:
                of = offset[8]
                am = amp[8]
                am2 = amp2[8]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]<110.:
                m1 = (offset[9]-offset[8])/(sim_freq[60]-sim_freq[50])
                b1 = offset[8]-m1*sim_freq[50]
                of = m1*sim_freq[f]+b1 
                m2 = (amp[9]-amp[8])/(sim_freq[60]-sim_freq[50])
                b2 = amp[8]-m2*sim_freq[50] 
                am = m2*sim_freq[f]+b2
                m3 = (amp2[9]-amp2[8])/(sim_freq[60]-sim_freq[50])
                b3 = amp2[8]-m3*sim_freq[50]
                am2 = m3*sim_freq[f]+b3                
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]==110.: 
                of = offset[9] 
                am = amp[9]
                am2 = amp2[9]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
    elif antenna=='100':
        for f in range(0,len(sim_freq)):
            if sim_freq[f]>80.:
                data[f] = 10**(sim_data[f]/10.)
            elif sim_freq[f]==80.:
                of = offset[6]
                am = amp[6]
                am2 = amp2[6]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]>70.:
                m1 = (offset[6]-offset[5])/(sim_freq[30]-sim_freq[20])
                b1 = offset[5]-m1*sim_freq[20]
                of = m1*sim_freq[f]+b1
                m2 = (amp[6]-amp[5])/(sim_freq[30]-sim_freq[20])
                b2 = amp[5]-m2*sim_freq[20]
                am = m2*sim_freq[f]+b2
                m3 = (amp2[6]-amp2[5])/(sim_freq[30]-sim_freq[20]) 
                b3 = amp2[5]-m3*sim_freq[20] 
                am2 = m3*sim_freq[f]+b3
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]==70.:
                of = offset[5]
                am = amp[5]
                am2 = amp2[5]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]>60.:
                m1 = (offset[5]-offset[4])/(sim_freq[20]-sim_freq[10])
                b1 = offset[4]-m1*sim_freq[10]
                of = m1*sim_freq[f]+b1
                m2 = (amp[5]-amp[4])/(sim_freq[20]-sim_freq[10])
                b2 = amp[4]-m2*sim_freq[10]
                am = m2*sim_freq[f]+b2
                m3 = (amp2[5]-amp2[4])/(sim_freq[20]-sim_freq[10]) 
                b3 = amp2[4]-m3*sim_freq[10] 
                am2 = m3*sim_freq[f]+b3
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]==60.:
                of = offset[4]
                am = amp[4]
                am2 = amp2[4]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]>50.:
                m1 = (offset[4]-offset[3])/(sim_freq[10]-sim_freq[0])
                b1 = offset[3]-m1*sim_freq[0]
                of = m1*sim_freq[f]+b1 
                m2 = (amp[4]-amp[3])/(sim_freq[10]-sim_freq[0])
                b2 = amp[3]-m2*sim_freq[0]
                am = m2*sim_freq[f]+b2
                m3 = (amp2[4]-amp2[3])/(sim_freq[10]-sim_freq[0]) 
                b3 = amp2[3]-m3*sim_freq[0] 
                am2 = m3*sim_freq[f]+b3
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))
            elif sim_freq[f]==50.:
                of = offset[3]
                am = amp[3]
                am2 = amp2[3]
                data[f] = (of+am*np.cos(az*2.))+am2*(1-np.cos(4.*az))

    return data

antenna = '70'
supdir = '../../supplemental_data/'
outdir = supdir

sim_data,sim_az,sim_alt,sim_freq = load_sim_files(antenna,supdir)
exp_data,exp_az,exp_alt = load_expand_files(antenna,supdir)
exp_freq = np.arange(20,120,10)
new_sd = 0.0000001*np.ones((len(sim_freq),len(sim_data[0])))
offset,amp = calc_sins(exp_data)
offset,amp,amp2 = calc_2nd_sin(antenna,sim_data,sim_alt,offset,amp)
for alt in range(0,len(sim_alt)):
    if sim_alt[alt]>=20*np.pi/180.:
            a = np.where(abs(sim_alt[alt]-exp_alt[0:91])<=0.01)[0]
            new_sd[:,alt] = lin_freq_fit(antenna,offset[:,a],amp[:,a],amp2[:,a],sim_freq,sim_data[:,alt],sim_az[alt])
#    else:
#        new_sd[:,alt] = 10**(sim_data[:,alt]/10.)

bad = np.where(new_sd<=0.)
new_sd[bad] = 0.0000001
new_sd = 10.*np.log10(new_sd)
np.save(supdir+'beam_simulations_ant70_50-110.npy',new_sd)

antenna = '100'
sim_data,sim_az,sim_alt,sim_freq = load_sim_files(antenna,supdir)
exp_data,exp_az,exp_alt = load_expand_files(antenna,supdir)
new_sd = 0.0000001*np.ones((len(sim_freq),len(sim_data[0])))
offset,amp = calc_sins(exp_data)
offset,amp,amp2 = calc_2nd_sin(antenna,sim_data,sim_alt,offset,amp)
for alt in range(0,len(sim_alt)):
    if sim_alt[alt]>=20*np.pi/180.:
        a = np.where(abs(sim_alt[alt]-exp_alt[0:91])<=0.01)[0]
        new_sd[:,alt] = lin_freq_fit(antenna,offset[:,a],amp[:,a],amp2[:,a],sim_freq,sim_data[:,alt],sim_az[alt])
#    else:
#        new_sd[:,alt] = 10**(sim_data[:,alt]/10.)

bad = np.where(new_sd<=0.)
new_sd[bad] = 0.0000001 
new_sd = 10.*np.log10(new_sd)
np.save(supdir+'beam_simulations_ant100_50-110.npy',new_sd)



