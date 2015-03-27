from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import ephem as eph
import numpy.polynomial.polynomial as poly


def find_gain(f,phi,theta,f_array,phi_array,theta_array,gaindb):
    """
    Determining the specific gain value to use from the beam data.
    f, phi, and theta assign the specific point of interest.
    f_array, phi_array,theta_array are the entire reference datasets
    gaindb is the gain array
    """
    
    f_index = where(f_array<=f)[0][-1]
    t_index = where(theta_array<=theta)[0][-1]
    p_index = where(phi_array<=phi)[0][-1]
    
    P = gaindb[t_index + MAXT(p_index+MAXP*f_index)]
    p = power(10.,0.1*P)
    
    return p


def gain_linear(t_index,p_index,f_index,gaindb):
    """
    Converts a gain value to linear scale.
    """
    P = gaindb[t_index+MAXT*(p_index+MAXP*f_index)]
    p = power(10.,0.05*P)

    return p

def find_fpt(f, phi,theta,f_array,phi_array,theta_array):
    """
    Determines the appropriate indices needed for gain data selection.
   
    Using find_fpt, then gain_linear is roughly the same as find_gain except for a square root in the power conversion. 
    """
    f_idx = where(f_array<=f)[0][-1]
    t_idx = where(theta_array<=theta)[0][-1]
    p_idx = where(phi_array<=phi)[0][-1]

    return f_idx,t_idx,p_idx

def findg(freq,phi,theta,gaindb,f,p,t):
    """
    calculates gain for the appropriate, f,p,t
    """
    epsilon = 1.e-3
    th1 = 90.-t
    
    if t>90.-epsilon:
        t -= epsilon
    if t < epsilon:
        t = epsilon
    
    if t > 90.:
        return 0.
    if t < 0.:
        return 0.

    if isnan(p)==True:
        p=180.
    elif p>180.:
        p -=180.
    elif p>90.:
        p = 180.-p
    elif p==0.:
        p += epsilon
    elif p==90.:
        p -= epsilon
   
    f_idx,p_idx,t_idx=find_fpt(f,p,t,freq,phi,theta)

   
    slope_t = (gain_linear(t_idx,p_idx-1,f_idx-1,gaindb)-gain_linear(t_idx-1,p_idx-1,f_idx-1,gaindb))/(theta[t_idx]-theta[t_idx-1])
    ans1 = gain_linear(t_idx-1,p_idx-1,f_idx-1,gaindb)+slope_t*(t-theta[t_idx-1])
  
    slope_t = (gain_linear(t_idx,p_idx,f_idx-1,gaindb)-gain_linear(t_idx-1,p_idx,f_idx-1,gaindb))/(theta[t_idx]-theta[t_idx-1])
    ans2 = gain_linear(t_idx-1,p_idx,f_idx-1,gaindb)+slope_t*(t-theta[t_idx-1])
  
    slope_p = (ans2-ans1)/(phi[p_idx]-phi[p_idx-1])
    ans3 = ans1+(slope_p*(p-phi[p_idx-1]))

    slope_t = (gain_linear(t_idx,p_idx-1,f_idx,gaindb)-gain_linear(t_idx-1,p_idx-1,f_idx,gaindb))/(theta[t_idx]-theta[t_idx-1])
    ans1 = gain_linear(t_idx-1,p_idx-1,f_idx,gaindb)+slope_t*(t-theta[t_idx-1]) 
   
    slope_t = (gain_linear(t_idx,p_idx,f_idx,gaindb)-gain_linear(t_idx-1,p_idx,f_idx,gaindb))/(theta[t_idx]-theta[t_idx-1])
    ans2 = gain_linear(t_idx-1,p_idx,f_idx,gaindb)+slope_t*(t-theta[t_idx-1])
   
    slope_p = (ans2-ans1)/(phi[p_idx]-phi[p_idx-1]) 
    ans4 = ans1+(slope_p*(p-phi[p_idx-1]))

    slope_f = (ans4-ans3)/(freq[f_idx]-freq[f_idx-1])
    ans5 = ans3 + slope_f*(f-freq[f_idx-1])

    return ans5

def galaxy_sky_map(direct,freq):
    """
    """
    file = loadtxt(direct+'radec'+str(int(freq))+'.dat')
    
    return file

def get_ra_dec(filename):
    """
    gets the appropriate ra,dec array for the angles.dat file
    (ra, dec for each healpix position). 
    """
    file = loadtxt(filename)
    ra_array = file[:,1]
    dec_array = 90.*ones(len(file))-file[:,0]

    return ra_array,dec_array

def antenna_beam_pattern(filename):
    """
    Converts the beam_simulations50-90.dat file into a gain array, with corresponding theta,phi and freq arrays.
    """
    
    freq = []
    theta = []
    phi = []
    gaindb = []
   
    f = open(filename)

    for line in range(0,21,1):
        singlef = f.readline()
#        freq_array.append(float(singlef.split(' ')[0]))
        for phi in range(0,181,1):
            for theta in range(0,181,1):
                singleline= f.readline()
                theta.append(float(singleline.split('  ')[0]))
                phi.append(float(singleline.split('  ')[1]))
                gaindb.append(float(singleline.split('  ')[2]))
                freq.append(float(singlef.split(' ')[0]))
    f.close()
    
    return freq,theta,phi,gaindb

def galaxy_temperature(lst,f,ras,decs,freqs,thetas,phis,gaindb):
    """
    """
    offset = -15.
    the_ras = ras*pi/180.
    the_decs = decs*pi/180.
    

    return
