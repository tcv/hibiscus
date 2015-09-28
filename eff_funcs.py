import numpy as np
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
    real = np.zeros(len(ant_data))
    imag = np.zeros(len(ant_data))
    freq = np.zeros(len(ant_data))
    for i in range(0,len(ant_data)):
        real[i] = [i][1]
        imag[i] = ant_data[i][2]
        freq[i] = ant_data[i][0]
    gamma = sqrt(real**2+imag**2)
    phi_gamma = arctan2(imag, real)
    delta_phi = cable_len*1.25E-10*2*np.pi*freq
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
    real_data = np.zeros(len(data.s))
    imag_data = np.zeros(len(data.s))
    freq_data = np.zeros(len(data.s))
    mag_data = np.zeros(len(data.s))
    phase_data = np.zeros(len(data.s))
    for i in range(0,len(data.s)):
        real_data[i] = np.real(data.s[i][0][0])
        imag_data[i] = np.imag(data.s[i][0][0])
        freq_data[i] = data.f[i]
        mag_data[i] = abs(data.s[i][0][0])
        phase_data[i] = np.angle(data.s[i][0][0])

    delta_phi = time_delay*2*np.pi*freq_data
    real_mod = mag_data*np.cos(phase_data+delta_phi)
    imag_mod = mag_data*np.sin(phase_data+delta_phi)

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
                E = np.sqrt(4*abs(Rant_func(F_ant)*Ramp_func(F_ant))/((Rant_func(F_ant)+Ramp_func(F_ant))**2+(Xant_func(F_ant)+Xamp_func(F_ant))**2))
            else:
                E = np.sqrt(4*abs(Rant_func(F_amp)*Ramp_func(F_amp))/((Rant_func(F_amp)+Ramp_func(F_amp))**2+(Xant_func(F_amp)+Xamp_func(F_amp))**2))
    elif len(R_ant)==1.0:
        if len(R_amp)>1.0:
            Ramp_func = itp.UnivariateSpline(F_amp,R_amp)
            Xamp_func = itp.UnivariateSpline(F_amp,X_amp)
            if F_amp[0]<=F_ant[0]:
                E = np.sqrt(4*abs(R_ant*Ramp_func(F_ant))/((R_ant+Ramp_func(F_ant))**2+(X_ant+Xamp_func(F_ant))**2))
            else:
                E = np.sqrt(4*abs(R_ant*Ramp_func(F_amp))/((R_ant+Ramp_func(F_amp))**2+(X_ant+Xamp_func(F_amp))**2))
        elif len(R_amp)==1.0:
            if F_amp[0]<=F_ant[0]:
                length = np.ones(len(F_ant))
                E = length*np.sqrt(4*abs(R_ant*R_amp)/((R_ant+R_amp)**2+(X_ant+X_amp)**2))
            else:
                length = np.ones(len(F_amp))
                E = length*np.sqrt(4*abs(R_ant*R_amp)/((R_ant+R_amp)**2+(X_ant+X_amp)**2))

    return E

