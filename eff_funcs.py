from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import ephem as eph

def imped(ant_data,cable_len):
    """
    Calculates impedence from the smith chart data.
    Assumes the data input is a three column array with
    real smith, imag smith, freq
    Can add a phase shift if needed.
    Also outputs the frequency in MHz.
    """
    real = []
    imag = []
    freq = []
    for i in range(0,len(ant_data)):
        real.append(ant_data[i][1])
        imag.append(ant_data[i][2])
        freq.append(ant_data[i][0])
    real = array(real)
    imag = array(imag)
    freq = array(freq)
    gamma = sqrt(real**2+imag**2)
    phi_gamma = arctan2(imag, real)
    delta_phi = cable_len*1.25E-10*2*pi*freq
    gamma_r = gamma*cos(phi_gamma+delta_phi)
    gamma_i = gamma*sin(phi_gamma+delta_phi)

    R = 50*(1-gamma**2)/((1-gamma_r)**2+gamma_i**2)
    X = 50*2*gamma_i/((1-gamma_r)**2+gamma_i**2)
    F = freq/1000000.0
    
    return R,X,F

def imped_skrf(file_name,time_delay):
    """
    Calculates impedence from s2p or s1p file
    Can add a phase shift if needed, but expects time in seconds.
    Also outputs the frequency in MHz.
    """
    data = rf.Network(file_name)
    real_data = []
    imag_data = []
    freq_data = []
    mag_data = []
    phase_data = []
    for i in range(0,len(data.s)):
        real_data.append(real(data.s[i][0][0]))
        imag_data.append(imag(data.s[i][0][0]))
        freq_data.append(data.f[i])
        mag_data.append(absolute(data.s[i][0][0]))
        phase_data.append(angle(data.s[i][0][0]))

    freq_data = array(freq_data)
    imag_data = array(imag_data)
    real_data = array(real_data)
    mag_data = array(mag_data)
    phase_data = array(phase_data)
    
    delta_phi = time_delay*2*pi*freq_data
    real_mod = mag_data*cos(phase_data+delta_phi)
    imag_mod = mag_data*sin(phase_data+delta_phi)

    R = 50*(1-mag_data**2)/((1-real_mod)**2+imag_mod**2)
    X = 100*imag_mod/((1-real_mod)**2+imag_mod**2)
    F = freq_data/1e6

    return R,X,F
    
def effic(R_ant, X_ant, F_ant, R_amp, X_amp, F_amp):
    """
    Calculates efficiency for a given set of amp, ant data.
    If freq range is different, freq range of Efficiency is smaller range.
    """
    if len(R_ant)>1.0:
        Rant_func = itp.UnivariateSpline(F_ant,R_ant)
        Xant_func = itp.UnivariateSpline(F_ant,X_ant)
        if len(R_amp)>1.0:
            Ramp_func = itp.UnivariateSpline(F_amp,R_amp)
            Xamp_func = itp.UnivariateSpline(F_amp,X_amp)
            if F_amp[0]<=F_ant[0]:
                E = sqrt(4*abs(Rant_func(F_ant)*Ramp_func(F_ant))/((Rant_func(F_ant)+Ramp_func(F_ant))**2+(Xant_func(F_ant)+Xamp_func(F_ant))**2))
            else:
                E = sqrt(4*abs(Rant_func(F_amp)*Ramp_func(F_amp))/((Rant_func(F_amp)+Ramp_func(F_amp))**2+(Xant_func(F_amp)+Xamp_func(F_amp))**2))
    elif len(R_ant)==1.0:
        if len(R_amp)>1.0:
            Ramp_func = itp.UnivariateSpline(F_amp,R_amp)
            Xamp_func = itp.UnivariateSpline(F_amp,X_amp)
            if F_amp[0]<=F_ant[0]:
                E = sqrt(4*abs(R_ant*Ramp_func(F_ant))/((R_ant+Ramp_func(F_ant))**2+(X_ant+Xamp_func(F_ant))**2))
            else:
                E = sqrt(4*abs(R_ant*Ramp_func(F_amp))/((R_ant+Ramp_func(F_amp))**2+(X_ant+Xamp_func(F_amp))**2))
        elif len(R_amp)==1.0:
            if F_amp[0]<=F_ant[0]:
                length = ones(len(F_ant))
                E = length*sqrt(4*abs(R_ant*R_amp)/((R_ant+R_amp)**2+(X_ant+X_amp)**2))
            else:
                length = ones(len(F_amp))
                E = length*sqrt(4*abs(R_ant*R_amp)/((R_ant+R_amp)**2+(X_ant+X_amp)**2))

    return E

