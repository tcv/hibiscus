from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import scipy.optimize as opt
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
        sub_data = zeros((len(time),len(freq)))
        sub_mask = zeros((len(time),len(freq)))
        num_mean = 0.
        for j in range(0,len(time)):
            if abs(stack_time[i]-time[j])<=(stack_time[1]-stack_time[0])/2.:
                sub_data[num_mean] = data[j]
                sub_mask[num_mean] = mask[j]
                num_mean = num_mean+1.
        sub_data = sub_data[0:num_mean]
        sub_mask = sub_mask[0:num_mean]
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
    lim_stack = zeros((len(data),len(data[0])))
    lim_mask = zeros((len(data),len(data[0])))
    lim_gsm = zeros((len(data),len(data[0])))
    lim_time = zeros(len(data))
    int = 0
    
    for i in range(0,len(data)):
        if sum(data[i])>0.:
            lim_stack[int] = data[i]
            lim_mask[int] =mask[i] 
            single_smooth = itp.UnivariateSpline(gsm_freq,gsm_data[i])
            lim_gsm[int] = single_smooth(freq)
            lim_time[int] =  gsm_time[i]
            int +=1
    lim_stack = lim_stack[0:int]
    lim_gsm = lim_gsm[0:int]
    lim_mask = lim_mask[0:int]
    lim_time = lim_time[0:int]

    return lim_stack, lim_mask, lim_gsm, lim_time

def time_mean(data,mask):
    """
    Calculates time mean for a 2 d array (1st ind time, 2nd ind freq).
    """
    mean_data = zeros(len(data[0]))
    mean_mask = zeros(len(data[0]))
    data = array(data)
    mask = array(mask)
    print shape(data),shape(mask)
    num_time = len(mask)
    num_freq = len(mask[0])
    mod_mask = zeros((num_time,num_freq))
    int = 0

    for i in range(0,num_time):
        single_mask = mask[i]
        if sum(single_mask)>=len(single_mask)/2.:
            new_mask = ones(len(single_mask))
        else:
            new_mask = zeros(len(single_mask))
        mod_mask[int] = new_mask
        int +=1

#    mod_mask = mod_mask[0:int]
    mod_mask = array(mod_mask)
    int2 = 0 
    for i in range(0,num_freq):
        single_mask = mod_mask[:,i]
        bad_num = sum(single_mask)
        if num_time<=bad_num:
            mean_data[int2] = 0.0
            mean_mask[int2] = 1.0
#            int2+=1
        else:
            single = ma.array(data[:,i],mask=single_mask)
            single_comp = ma.compressed(single)
            mean_data[int2]= ma.mean(single_comp)
            mean_mask[int2] = 0.0
        int2+=1
    
#    mean_data = mean_data[0:int2]
#    mean_mask = mean_mask[0:int2]

    return mean_data,mean_mask

def gain_calc(data,masked,gsm,K0):
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

def poly_fore(data,masked,freq,minf,maxf,n,std):
    """
    Calculates a polynomial fit for data.
    Inputs:
    data - single frequency dependent spectrum
    masked - corresponding mask
    freq - corresponding frequency array
    minf, maxf - min and max freq if want to truncate range 
    n - index of polynomial fitting
    std - 1/weights, set to ones if not needed.
    Output:
    dfit - polynomial fit spectrum
    fit_params - parameters for the fit
    """
    
    data_array = ma.array(data,mask=masked)
    data_comp = ma.compressed(data_array)
    freq_array = ma.array(freq,mask=masked)
    freq_comp = ma.compressed(freq_array)
    std_comp = ma.compressed(ma.array(std, mask=masked))


    min_ind = 0
    max_ind = -1
    if minf>freq_comp[0]:
        min_ind = where(freq_comp<=minf)[0][-1]
    if maxf<freq_comp[-1]:
        max_ind = where(freq_comp<=maxf)[0][-1]

    log_data = log10(data_comp[min_ind:max_ind])
    log_freq = log10(freq_comp[min_ind:max_ind]/70.)
    weights = 1/std_comp[min_ind:max_ind]
    
    fit_params = poly.polyfit(log_freq,log_data,n,w=weights)
    dfit = 10**(poly.polyval(log10(freq/70.),fit_params))

    return dfit, fit_params
