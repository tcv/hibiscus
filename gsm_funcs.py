"""
Set of scripts for doing things related to the gsm data generation process.
"""

from numpy import *
import pylab
import numpy
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import ephem as eph
import numpy.polynomial.polynomial as poly
from matplotlib import cm
import matplotlib.pyplot as plt

def get_gsm_radec(filename):
    """
    gets the appropriate ra,dec array for the angles.dat file
    (ra, dec for each healpix position).
    goes with the gsm data. 
    """
    file = loadtxt(filename)
    ra_array = file[:,1]
    dec_array = 90.*ones(len(file))-file[:,0]
    file = []

    ra_array = array(ra_array)
    dec_array = array(dec_array)

    return ra_array,dec_array

def antenna_beam_pattern(filename,input_freqs):
    """
    Converts the beam_simulations50-90.dat file into a gain array,
    with corresponding theta and phi arrays.
    Done for a given frequency
    """
    
    theta_array = zeros(181*181)
    phi_array = zeros(181*181)
    gaindb = -50.*ones(181*181)
    gaindb2 = -50.*ones(181*181)
    larger = 0
    freq1 = 0
    freq2 = 0
 
    f = open(filename)

    for line in range(0,21,1):
        singlef = f.readline()
        for phi in range(0,181,1):
            for theta in range(0,181,1):
                singleline= f.readline()
                if (float(singlef.split('  ')[0])-float(input_freqs))<=0.0:
                    theta_array[phi*181+theta] = (float(singleline.split('  ')[0]))
                    phi_array[phi*181+theta] = (float(singleline.split('  ')[1]))
                    gaindb[phi*181+theta] = (float(singleline.split('  ')[2]))
                    freq1 = float(singlef.split('  ')[0])
                if (float(singlef.split('  ')[0])-float(input_freqs))==0.:
                    theta_array[phi*181+theta] = (float(singleline.split('  ')[0]))
                    phi_array[phi*181+theta] = (float(singleline.split('  ')[1]))
                    gaindb[phi*181+theta] = (float(singleline.split('  ')[2]))
                    freq1 = float(singlef.split('  ')[0])
                    freq2 = freq1
                if (float(singlef.split('  ')[0])-float(input_freqs))>0.:
                    if freq1!=freq2:
                        if larger<(181*181):
                            gaindb2[phi*181+theta] = (float(singleline.split('  ')[2]))
                            larger = larger + 1
                            freq2 = float(singlef.split('  ')[0]) 
    f.close()

    if freq2==0:
        freq2 = 130.

    if freq1==freq2:
        gain_act = gaindb
    else:
        m = (gaindb-gaindb2)/(freq1-freq2)
        b = (gaindb)-m*freq1
        gain_act = m*float(input_freqs)+b

#    print "For the Antenna Beam Simulation"
#    print "Input Frequency is: ",input_freqs, " MHz"
#    print "Used Frequencies are: ", freq1," and ", freq2, " MHz"

    return theta_array,phi_array,gain_act

def antenna_beam_full(filename,input_freqs):
    """
    expanding the simulated antenna beam
     to include data for altitude less than zero smoothly.
    """
    theta_sim,phi_sim,gaindb=antenna_beam_pattern(filename,input_freqs)
    
    theta_sim = theta_sim*pi/180.
    phi_sim = phi_sim*pi/180.
#    print 'Theta Range: ',ma.min(theta_sim),ma.max(theta_sim)
#    print 'Phi Range: ',ma.min(phi_sim),ma.max(phi_sim)
    curr_alt = zeros(len(phi_sim)*2)
    curr_az = zeros(len(phi_sim)*2)
    az_inds = zeros(181)
    azs = arange(0,2.*pi,2*pi/181.) 
    beam_data = -50.*ones(len(phi_sim)*2)
    for p in range(0,len(phi_sim)):
        if theta_sim[p]<=0:
            curr_az[p] = 2*pi-phi_sim[p]
            curr_alt[p] = pi/2.+theta_sim[p]
            curr_az[-p] = curr_az[p]
            curr_alt[-p] = -1.*curr_alt[p]
        else:
            curr_az[p] = phi_sim[p]
            curr_alt[p] = pi/2.-theta_sim[p]
            curr_az[-p] = curr_az[p] 
            curr_alt[-p] = -1.*curr_alt[p]
        if curr_alt[p]==0:
            sin_ind = where(abs(curr_az[p]-azs)<(2*pi/181.))
            az_inds[sin_ind] = p
      
        beam_data[p] = gaindb[p]

    for p in range(len(phi_sim),2*len(phi_sim)):
        for th in range(0,181):
            if abs(curr_az[p]-azs[th])<=(pi/180.):
                beam_data[p] = gaindb[az_inds[th]]
        
    return curr_az,curr_alt,beam_data

def gsm_temps(gsmdir,input_freqs,inds):
    """
    loads the gsm data into a temperature array,
    with corresponding freq array (ra/dec arrays from get_gsm_radec). 
    """
    
    if float(input_freqs)>=50.:
        if float(input_freqs)<=112.:
            if input_freqs%1.<=0.01:    
                freqs1 = int(input_freqs)
                fname = gsmdir+'radec'+str(int(freqs1))+'.dat'
                gsmdata = loadtxt(fname)
                freqs2 = freqs1
            else:
                freqs1 = float(int(input_freqs))
                freqs2 = freqs1+1
                fname1 = gsmdir+'radec'+str(int(freqs1))+'.dat'
                gsmdata1 = loadtxt(fname1)
                fname2 = gsmdir+'radec'+str(int(freqs2))+'.dat'
                gsmdata2 = loadtxt(fname2)
                m = (gsmdata1-gsmdata2)/(freqs1-freqs2)
                b = gsmdata1-m*freqs1
                gsmdata = m*float(input_freqs)+b
        else:
            freqs1 = 130.
            freqs2 = freqs1
            fname = gsmdir+'radec112.dat'
            gsmdatat = loadtxt(fname)
            gsmdata = zeros(len(gsmdatat))
    else:
        freqs1 = 40.
        freqs2 = freqs1
        fname = gsmdir+'radec50.dat'
        gsmdatat = loadtxt(fname)
        gsmdata = zeros(len(gsmdatat)) 
    
#    print "For the GSM Model Data"
#    print "Input Frequency is: ", input_freqs, " MHz"
#    print "Used Frequencies are: ",freqs1, " and " ,freqs2, " MHz"
    used_data = zeros(len(inds))
    for i in range(0,len(inds)):
        used_data[i] =gsmdata[inds[i]]

    return used_data

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
    point = eph.FixedBody()
    point._ra = curr_ra
    point._dec = curr_dec
    point.compute(site)
    cur_alt = point.alt
    cur_az = point.az
    return cur_alt, cur_az

def radec_sim(curr_az,curr_alt,lat,lon,elevation,time,idate):
    """
    Calculate ra,dec using lst, lat, az and alt
    ra runs from 0 to 2 pi
    dec runs from -pi/2 to pi/2
    """

    site = eph.Observer()
    site.lon = lon
    site.lat = lat
    site.elevation = elevation
    date = eph.date(idate)+time/24.
    site.date = date
    site.pressure =0

#    print 'Altitude range: ',ma.min(curr_alt), ma.max(curr_alt)
#    print 'Azimuth range: ',ma.min(curr_az),ma.max(curr_az)
     
    ra = zeros(len(curr_az))
    dec = zeros(len(curr_az))
    for p in range(0,len(curr_az)):
        ra[p],dec[p]=site.radec_of(curr_az[p],curr_alt[p])

#    print 'RA range is: ',ma.min(ra),ma.max(ra)
#    print 'DEC range is: ',ma.min(dec),ma.max(dec)

    sim_var = vstack((ra,dec)).T
    
    return sim_var
 
def calc_ra_shift(lat,lon,elevation,time1,time2,idate):
    """
    Calculate the amount of shift in RA needed between times.
    """
    site = eph.Observer()
    site.lon = lon
    site.lat = lat
    site.elevation = elevation
    date1 = eph.date(idate)+time1/24.
    date2 = eph.date(idate)+time2/24.
    site.pressure =0

    site.date = date1
    ra1,dec1 = site.radec_of(pi/2.,pi/4.)
    site.date = date2
    ra2,dec2 = site.radec_of(pi/2.,pi/4.)
    radiff = float(int((ra2-ra1)*180./pi*1e2)/1e2)
    return radiff

def gsm_comp(gsmdata,ras,decs,lat,lon,elevation,time,idate):
    """
    Creates the appropriate gsm array for a given time and freq
    """

    alts = zeros(len(ras))
    azs = zeros(len(ras))
    #freqs_gsm = []
    gsm_array = zeros(len(ras))
    len_a  = 0
    for i in range(0,len(ras)):
        sin_alt, sin_az = azel_loc(ras[i],decs[i],lat,lon,elevation,time,idate)
        if sin_alt>=0:
            alts[len_a] = sin_alt
            azs[len_a] = sin_az
            gsm_array[len_a] = gsmdata[i]
            len_a = len_a+1

        
    alts = array(alts[0:len_a-1])
    azs = array(azs[0:len_a-1])
    gsm_array = array(gsm_array[0:len_a-1])

    gsm_var = vstack((alts,azs)).T
    alts = []
    azs = []    
    
    return gsm_array, gsm_var


def sim_comp(beamfile,freqs):
    """
    Generates the sim dataset in the right shape for a given freq.
    """

    theta_sim,phi_sim,gaindb = antenna_beam_pattern(beamfile,freqs)

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

def sim_beam_interp(gsm_var,gaindb,sim_var):
    """
    Generates the interpolated sim beam for the GSM ra/decs
    """
# Alternate interpolator, but takes longer to run
#    func_gain = itp.CloughTocher2DInterpolator(sim_var*180./pi,gaindb)

#Simplest linear interpolator. 
    func_gain = itp.LinearNDInterpolator(sim_var*180./pi,gaindb)


    gsm_gain = func_gain(gsm_var)
    gain_beam = pow(10.,0.05*gsm_gain)
    nandata = where(isnan(gain_beam))
    for i in range(0,len(nandata[0])):
        gain_beam[nandata[0][i]]=0.0

    return gain_beam

def ant_beam_simple(gsm_array,gain_beam):
    """
    Just doing the calculation of the single beam value for two arrays.
    Both arrays must have the same RA/DEC array. 
    """    
    full_beam = gain_beam*gsm_array
    nandata = where(isnan(full_beam))
    for i in range(0,len(nandata[0])):
        full_beam[nandata[0][i]]=0.0
 
    summed_beam = ma.sum(ma.sum(full_beam,axis=0),axis=0)
    summed_sim = ma.sum(ma.sum(gain_beam,axis=0),axis=0)
    final_result = summed_beam/summed_sim

    return final_result

def ant_beam(gsm_array, gsm_var,gaindb,sim_var,label,freq,plotf):
    """
    Combines the gsm and sim datasets for a given place/time.
    Note I've limited the frequency range that is loaded to avoid memory errors
    Re-wrote to limit to a single frequency 
    Expects sim_var,gsm_var to be in degrees. 
    """
    
    gain_beam = sim_beam_interp(gsm_var,gaindb,sim_var)
    full_beam = gain_beam*gsm_array

    nandata = where(isnan(full_beam))
    for i in range(0,len(nandata[0])):
        full_beam[nandata[0][i]]=0.0

    summed_beam = ma.sum(ma.sum(full_beam,axis=0),axis=0)
    summed_sim = ma.sum(ma.sum(gain_beam,axis=0),axis=0)

#Allows you to make plots to check results at a single frequency only if you set plotf to be within the frequency range of the data.     
    if freq==plotf:
        plt.rc('font',size=8)
        plt.subplot(411)
        plt.scatter(sim_var[:,0]*180./pi,sim_var[:,1]*180./pi,s=1,linewidth=0,c=pow(10.,0.05*gaindb),vmin=0,vmax=3,cmap=cm.jet)
        plt.colorbar() 
        plt.xlim(0,360)
        plt.ylim(-90,90)
        plt.ylabel('DEC (degrees)')
        plt.title('Simulated HIbiscus Beam (linear power)')

        plt.subplot(412)
        plt.scatter(gsm_var[:,0],gsm_var[:,1],s=1,linewidth=0,c=gain_beam,vmin=0,vmax=3,cmap=cm.jet)
        plt.colorbar()
        plt.xlim(0,360)
        plt.ylim(-90,90)
        plt.ylabel('DEC (degrees)')
        plt.title('Interpolated HIbiscus Beam (linear power)')

        plt.subplot(413)
        plt.scatter(gsm_var[:,0],gsm_var[:,1],s=1,linewidth=0,c=gsm_array,vmin=0,vmax=2e4,cmap=cm.jet)
        plt.colorbar()
        plt.xlim(0,360)
        plt.ylim(-90,90)
        plt.ylabel('DEC (degrees)')
        plt.title('GSM Data (Kelvin)')

        plt.subplot(414) 
        plt.scatter(gsm_var[:,0],gsm_var[:,1],s=1,linewidth=0,c=full_beam,vmin=0,vmax=5e4,cmap=cm.jet)
        plt.colorbar() 
        plt.xlim(0,360)
        plt.ylim(-90,90)
        plt.xlabel('RA (degrees)')  
        plt.ylabel('DEC (degrees)')
        plt.title('Expected Signal (Kelvin)')
        plt.subplots_adjust(hspace=0.4)
        plt.savefig(label,dpi=300)
        plt.clf()

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


