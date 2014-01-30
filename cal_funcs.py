from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import ephem as eph
import numpy.polynomial.polynomial as poly

def match_binning(gsm_time,freq,time,data,mask):
    """
    Process for matching the time binning of two datasets.
    Initially designed to match signal data to gsm datasets.
    Assumes that GSM data has coarser resolution.
    Inputs:
    gsm_time - time that is being matched to.
    freq - frequency array (only need length)
    time - original time array of data
    data - original data array
    mask - original mask array

    outputs are new data and mask array
    """
    stack_time = gsm_time
    stack_data = zeros((len(stack_time),len(freq)))
    stack_mask = zeros((len(stack_time),len(freq)))
    for i in range(0,len(stack_time)):
        sub_data = []
        sub_mask = []
        num_mean = 0.
        for j in range(0,len(time)):
            if abs(stack_time[i]-time[j])<=(stack_time[1]-stack_time[0])/2.:
                sub_data.append(data[j])
                sub_mask.append(mask[j])
                num_mean = num_mean+1.
        sub_data = array(sub_data)
        sub_mask = array(sub_mask)
        if num_mean>=1.0:
            for f in range(0,len(freq)):
                if sum(sub_mask[:,f])==len(sub_mask[:,f]):
                    stack_data[i,f]=0.0
                    stack_mask[i,f]=1.0
                else:
                    single_data = ma.array(sub_data[:,f],mask=sub_mask[:,f])
                    single_comp = ma.compressed(single_data)
                    stack_data[i,f] = ma.mean(single_comp)
    

    return stack_data, stack_mask

def lim_bin(freq,data,mask,gsm_freq,gsm_data,gsm_time):
    """
    A module to separate a set of matching arrays into only 
    the portion of the array that has data for both arrays.
    Necessary for data when we don't have full days.
    Inputs:
    freq - frequency array for data
    data - data array
    mask - mask array for data
    gsm_freq - frequency array for gsm data
    gsm_data - gsm data array
    gsm_time - gsm time array
    
    Outputs:
    lim_stack = new data array
    lim_mask = matching mask array
    lim_gsm = matching gsm array
    lim_time = matching time array
    """
    lim_stack = []
    lim_mask = []
    lim_gsm = []
    lim_time = []
    
    for i in range(0,len(data)):
        if sum(stack_data[i])>0.:
            lim_stack.append(data[i])
            lim_mask.append(mask[i])
            single_smooth = itp.UnivariateSpline(gsm_freq,gsm_data[i])
            lim_gsm.append(single_smooth(freq))
            lim_time.append(gsm_time[i])
    lim_stack = array(lim_stack)
    lim_gsm = array(lim_gsm)
    lim_mask = array(lim_mask)
    lim_time = array(lim_time)

    return lim_stack, lim_mask, lim_gsm, lim_time

def time_mean(data,mask):
    """
    Calculates time mean for a 2 d array (1st ind time, 2nd ind freq).
    """
    mean_data = []
    mean_mask = []
    for i in range(0,len(data[0])):
        if len(mask[:,i])==sum(mask[:,i]):
            mean_data.append(0.0)
            mean_mask.append(1.0)
        else:
            single = ma.array(data[:,i],mask=mask[:,i])
            single_comp = ma.compressed(single)
            mean_data.append(ma.mean(single_comp))
            mean_mask.append(0.0)
    
    mean_data = array(mean_data)
    mean_mask = array(mean_mask)

    return mean_data,mean_mask

def gain_calc(data,masked,gsm,K0)
    """
    Calculates the gain using least squares for a single frequency.
    Inputs:
    data - single freq array of data
    gsm  - similar single freq array of gsm data
    masked - mask for data
    K0 - preliminary guess for gain
    """
    fK = lambda K,d,g: K*d-g

    d0_array = ma.array(data,mask=masked)
    d0 = ma.compressed(d0_array)
    g0_array = ma.array(gsm,mask=masked)
    g0 = ma.compressed(g0_array)
    Kf = opt.leastsq(fK,K0,args=(d0,g0),maxfev=100000)
    
    return Kf[0]

def poly_fore(data,masked,freq,minf,maxf,n)
    """
    Calculates a polynomial fit for data.
    Inputs:
    data - single frequency dependent spectrum
    masked - corresponding mask
    freq - corresponding frequency array
    minf, maxf - min and max freq if want to truncate range 
    n - index of polynomial fitting
    Output:
    dfit - polynomial fit spectrum
    fit_params - parameters for the fit
    """
    
    data_array = ma.array(data,mask=masked)
    data_comp = ma.compressed(data_array)
    freq_array = ma.array(freq,mask=masked)
    freq_comp = ma.compressed(freq_array)

    min_ind = 0
    max_ind = -1
    if minf>freq_comp[0]:
        min_ind = where(freq_comp<=minf)[0][-1]
    if maxf<freq_comp[-1]:
        max_ind = where(freq_comp<=maxf)[0][-1]

    log_data = log10(data_comp[min_ind:max_ind])
    log_freq = log10(freq_comp[min_ind:max_ind]/70.)
    
    fit_params = poly.polyfit(log_freq,log_data,n)
    dfit = 10**(poly.polyval(log10(freq/70.),fit_params))

    return dfit, fit_params
