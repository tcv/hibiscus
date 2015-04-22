from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import ephem as eph
import numpy.polynomial.polynomial as poly

def get_gsm_radec(filename):
    """
    gets the appropriate ra,dec array for the angles.dat file
    (ra, dec for each healpix position).
    goes with the gsm data. 
    """
    file = loadtxt(filename)
    ra_array = file[:,1]
    dec_array = 90.*ones(len(file))-file[:,0]

    ra_array = array(ra_array)
    dec_array = array(dec_array)

    return ra_array,dec_array

def antenna_beam_pattern(filename,input_freqs):
    """
    Converts the beam_simulations50-90.dat file into a gain array,
    with corresponding theta,phi and freq arrays.
    """
    
    freq_array = []
    theta_array = []
    phi_array = []
    gaindb = []
   
    f = open(filename)

    for line in range(0,21,1):
        singlef = f.readline()
        for phi in range(0,181,1):
            for theta in range(0,181,1):
                singleline= f.readline()
                if (float(singlef.split('  ')[0])>=input_freqs):
                    if len(freq_array)<=181*181:
#                    if (float(singlef.split('  ')[0])<input_freqs+3.):
                        theta_array.append(float(singleline.split('  ')[0]))
                        phi_array.append(float(singleline.split('  ')[1]))
                        gaindb.append(float(singleline.split('  ')[2]))
                        freq_array.append(float(singlef.split(' ')[0]))
    f.close()

    freq_array = array(freq_array)
    theta_array = array(theta_array)
    phi_array = array(phi_array)
    gaindb = array(gaindb)
    
    return freq_array,theta_array,phi_array,gaindb

def gsm_temps(gsmdir,input_freqs):
    """
    loads the gsm data into a temperature array,
    with corresponding freq array (ra/dec arrays from get_gsm_radec). 
    """

    
    freqs = int(input_freqs)

    gsmdata = []
    fname = gsmdir+'radec'+str(int(freqs))+'.dat'
    single = loadtxt(fname)
    gsmdata.append(single)

    gsmdata = array(gsmdata)

    return freqs, gsmdata


def azel_loc(ra,dec,lat,lon,elevation,time,idate):
    """
    Uses ephem to convert ra/dec to az/el for a given location.
    alt runs from -pi/2 to pi/2
    az  runs from 0 to 2 pi
    """
    
    site = eph.Observer()
    site.lon = lon
    site.lat = lat
    site.elevation = elevation
    date = eph.date(idate)+time/24.
    site.date = date
    site.pressure =0

    curr_ra = eph.degrees(ra*pi/180.)
    curr_dec = eph.degrees(dec*pi/180.)
    db_entry = "curr_pt,f/J,"+str(curr_ra)+","+str(curr_dec)+",-1"
    point = eph.readdb(db_entry)

    point.compute(site)

    cur_alt = point.alt
    cur_az = point.az

    return cur_alt, cur_az

# At this point I have two data arrays, for the simulation and gsm data.
# Now I need code to calculate a expected temperature given
# Frequency and time.

def gsm_comp(gsmdata,ras,decs,lat,lon,elevation,time,idate):
    """
    Creates the appropriate gsm array for a given time and freq
    """

    alts = []
    azs = []
    #freqs_gsm = []
    gsm_array = []
    for i in range(0,len(ras)):
        sin_alt, sin_az = azel_loc(ras[i],decs[i],lat,lon,elevation,time,idate)
        if sin_alt>=0:
            alts.append(sin_alt)
            azs.append(sin_az)
    #        freqs_gsm.append(freq_gsm)
            gsm_array.append(gsmdata[0,i])
            
    alts = array(alts)
    azs = array(azs)
    #freqs_gsm = array(freqs_gsm)
    gsm_array = array(gsm_array)

    gsm_var = vstack((alts,azs)).T
    
    return gsm_array, gsm_var

def sim_comp(beamfile,freqs):
    """
    Generates the sim dataset in the right shape for a given freq.
    """

    freqs_sim,theta_sim,phi_sim,gaindb = antenna_beam_pattern(beamfile,freqs)

    theta_sim = theta_sim*pi/180.
    phi_sim = phi_sim*pi/180.

    alt_sim = zeros(len(phi_sim))
    az_sim = zeros(len(phi_sim))
    for p in range(0,len(phi_sim)):
        if theta_sim[p]<0:
            az_sim[p] = 2*pi - phi_sim[p]
            alt_sim[p] = pi/2. +theta_sim[p]
        else:
            az_sim[p] = phi_sim[p]
            alt_sim[p] = pi/2.-theta_sim[p]
            
    sim_var = vstack((alt_sim,az_sim)).T

    return gaindb, sim_var
    

def ant_beam(gsm_array, gsm_var, gaindb, sim_var):
    """
    Combines the gsm and sim datasets for a given place/time.
    Note I've limited the frequency range that is loaded to avoid memory errors
    Re-wrote to limit to a single frequency 
    """

    grid_alt, grid_az = mgrid[0:pi/2.:90j,0:2.*pi:180j]

    grid_gain = itp.griddata(sim_var,gaindb,(grid_alt,grid_az),method='linear')
    grid_temp = itp.griddata(gsm_var,gsm_array,(grid_alt,grid_az),method='linear')

    full_beam = pow(10.,0.05*grid_gain)*grid_temp
    
    summed_beam = ma.sum(ma.sum(full_beam,axis=0),axis=0)
    summed_sim = ma.sum(ma.sum(pow(10.,0.05*grid_gain),axis=0),axis=0)
#    print shape(summed_beam)

    final_result = summed_beam/summed_sim

    return final_result
    

def find_fpt(f, phi,theta,f_array,phi_array,theta_array):
    """
    Determines the appropriate indices needed for beam sim data selection.
    """
    f_idx = where(f_array<=f)[0][-1]
    t_idx = where(theta_array<=theta)[0][-1]
    p_idx = where(phi_array<=phi)[0][-1]

    return f_idx,t_idx,p_idx

def gain_linear(t_index,p_index,f_index,gaindb):
    """
    Converts a gain value to linear scale.
    Question of factor of 1/2 in power for sim data.
    """

    MAXT = 181.
    MAXP = 181.

    P = gaindb[t_index+MAXT*(p_index+MAXP*f_index)]
    p = power(10.,0.05*P)

    return p

def findg(freq,phi,theta,gaindb,f,p,t):
    """
    Interpolates the beam sim data for any f, p, t indices from find_fpt
    as long as as the values are within the ranges of the 
    freq, theta and phi arrays. 
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

beamfile = '../beam_simulations50-90.dat'
radecfile = '../angles.dat'
gsmdir = '../galaxy_maps_radec/'
lat = '-30.727206'
lon = '21.479055'
elevation = 1080
time = 3.5*24.
idate = '2015/04/01'
freqs = 60
