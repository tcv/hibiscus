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
#    dec_array = file[:,0]
    file = []

    ra_array = array(ra_array)
    dec_array = array(dec_array)

    return ra_array,dec_array

def antenna_beam_pattern(filename,input_freqs):
    """
    Converts the beam_simulations50-90.dat file into a gain array,
    with corresponding theta,phi and freq arrays.
    """
    
    freq_array = zeros(181*181)
    theta_array = zeros(181*181)
    phi_array = zeros(181*181)
    gaindb = -50.*ones(181*181)
    theta_array2 = zeros(181*181)
    phi_array2 = zeros(181*181)
    gaindb2 = -50.*ones(181*181)
    larger = 0
    freq1 = 0
    freq2 = 0

#    print len(theta_array), 181*181
   
    f = open(filename)

    for line in range(0,21,1):
        singlef = f.readline()
#        print float(singlef.split('  ')[0]),float(input_freqs)
        for phi in range(0,181,1):
            for theta in range(0,181,1):
                singleline= f.readline()
#                if (float(singlef.split('  ')[0])>=float(input_freqs)):
#                if abs(float(singlef.split(' ')[0])-float(input_freqs))<1.0:
#                    if (float(singlef.split('  ')[0])<input_freqs+3.):
                if (float(singlef.split('  ')[0])-float(input_freqs))<=0.0:
                    theta_array[phi*181+theta] = (float(singleline.split('  ')[0]))
                    phi_array[phi*181+theta] = (float(singleline.split('  ')[1]))
                    gaindb[phi*181+theta] = (float(singleline.split('  ')[2]))
                    freq1 = float(singlef.split('  ')[0])
#                    freq2 = freq1
#                    freq_array[phi*181+theta] = (float(singlef.split(' ')[0]))
                if (float(singlef.split('  ')[0])-float(input_freqs))==0.:
                    theta_array[phi*181+theta] = (float(singleline.split('  ')[0]))
                    phi_array[phi*181+theta] = (float(singleline.split('  ')[1]))
                    gaindb[phi*181+theta] = (float(singleline.split('  ')[2]))
                    freq1 = float(singlef.split('  ')[0])
                    freq2 = freq1
                if (float(singlef.split('  ')[0])-float(input_freqs))>0.:
                    if freq1!=freq2:
                        if larger<(181*181):
                            theta_array2[phi*181+theta] = (float(singleline.split('  ')[0]))
                            phi_array2[phi*181+theta] = (float(singleline.split('  ')[1]))
                            gaindb2[phi*181+theta] = (float(singleline.split('  ')[2]))
                            larger = larger + 1
                            freq2 = float(singlef.split('  ')[0]) 
#                    else:
#                        print larger
#                    freq_array[phi*181+theta] = (float(singlef.split(' ')[0]))

    f.close()
#    print sum(gaindb),sum(gaindb2)

#    if larger==0:
    if freq2==0:
        freq2 = 130.

    if freq1==freq2:
        gain_act = gaindb
    else:
        m = (gaindb-gaindb2)/(freq1-freq2)
        b = (gaindb)-m*freq1
        gain_act = m*float(input_freqs)+b


#    freq_array = array(freq_array)
#    theta_array = array(theta_array)
#    phi_array = array(phi_array)
#    gaindb = array(gaindb)


#    print ma.sum(gaindb)    
    print "For the Antenna Beam Simulation"
    print "Input Frequency is: ",input_freqs, " MHz"
    print "Used Frequencies are: ", freq1," and ", freq2, " MHz"

    return freq_array,theta_array,phi_array,gain_act

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
    
    print "For the GSM Model Data"
    print "Input Frequency is: ", input_freqs, " MHz"
    print "Used Frequencies are: ",freqs1, " and " ,freqs2, " MHz"
    used_data = zeros(len(inds))
    for i in range(0,len(inds)):
        used_data[i] =gsmdata[inds[i]]

    return freqs,used_data

def trunc_gsm(ras,decs,ra_var,dec_var):
    """
    Shrinks down the size of the gsmdata array
    """
    new_ras = zeros(ra_var*dec_var)
    new_decs = zeros(ra_var*dec_var)
    diff = (360.+180.)*ones(len(new_ras))
    new_gsm = zeros(ra_var*dec_var)
    rwid = ra_var/180.
    dwid = dec_var/360.
    rlist = arange(0,180+rwid,1/rwid)
    dlist = arange(0,360+dwid,1/dwid)
    for i in range(0,ra_var):
        for j in range(0,dec_var):
            new_decs[i*dec_var+j] = rlist[i]
            new_ras[i*dec_var+j] = dlist[j]
     
    for a in range(0,len(ras)):
        value = abs(new_ras-ras[a]*ones(len(new_ras)))+abs(new_decs-decs[a]*ones(len(new_decs)))
        ind = argmin(value)
        val = amin(value)
        if val<diff[ind]:
            diff[ind] = val
            new_gsm[ind] = a
   

    return new_ras,new_decs,new_gsm

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
#    db_entry = "curr_pt,f/J,"+str(curr_ra)+","+str(curr_dec)+",-1"
#    point = eph.readdb(db_entry)
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

    alts = zeros(len(ras))
    azs = zeros(len(ras))
    #freqs_gsm = []
    gsm_array = zeros(len(ras))
    len_a  = 0
    for i in range(0,len(ras)):
#    i = 1000
#    while 1:
        sin_alt, sin_az = azel_loc(ras[i],decs[i],lat,lon,elevation,time,idate)
        if sin_alt>=0:
            alts[len_a] = sin_alt
            azs[len_a] = sin_az
    #        freqs_gsm.append(freq_gsm)
            gsm_array[len_a] = gsmdata[i]
            len_a = len_a+1

        
    alts = array(alts[0:len_a-1])
    azs = array(azs[0:len_a-1])
    #freqs_gsm = array(freqs_gsm)
    gsm_array = array(gsm_array[0:len_a-1])

    gsm_var = vstack((alts,azs)).T
    alts = []
    azs = []    
    
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
    

def ant_beam(gsm_array, gsm_var, gaindb, sim_var,label):
    """
    Combines the gsm and sim datasets for a given place/time.
    Note I've limited the frequency range that is loaded to avoid memory errors
    Re-wrote to limit to a single frequency 
    """
    adj = 0.001
    grid_alt, grid_az = numpy.mgrid[0+adj:pi/2.-adj:180j,0+adj:2.*pi-adj:720j]

#    grid_gain = itp.griddata(sim_var,gaindb,(grid_alt,grid_az),method='cubic')
#    grid_temp = itp.griddata(gsm_var,gsm_array,(grid_alt,grid_az),method='cubic')
#    print amax(gsm_var[:,0]),amin(gsm_var[:,0])
#    print amax(gsm_var[:,1]),amin(gsm_var[:,0])
#    print shape(gsm_var)
    func_gain = itp.CloughTocher2DInterpolator(sim_var,gaindb)
    gsm_gain = func_gain(gsm_var)
    print shape(gsm_gain)
#    func_temp = itp.Rbf(gsm_var[:,0],gsm_var[:,1],gsm_array)
#    grid_temp = func_temp(grid_alt[:,0],grid_az[0,:])
#    grid_gain = func_gain(grid_alt[:,0],grid_az[0,:])
      
#    gain_beam = pow(10.,0.05*grid_gain)
#    full_beam = gain_beam*grid_temp
    gain_beam = pow(10.,0.05*gsm_gain)
    full_beam = gain_beam*gsm_array

    nandata = where(isnan(full_beam))
    for i in range(0,len(nandata[0])):
        full_beam[nandata[0][i]]=0.0
        gain_beam[nandata[0][i]]=0.0

    print shape(where(isnan(full_beam)))

    summed_beam = ma.sum(ma.sum(full_beam,axis=0),axis=0)
    summed_sim = ma.sum(ma.sum(gain_beam,axis=0),axis=0)
    print amax(full_beam/summed_sim)
    

#    n = plt.Normalize(0.,90.)
#    f,(plt1,plt2,plt3) = plt.subplots(3,1,sharex=True)
#    f.xlim(0,360)
    plt.subplot(311)
#    pylab.imshow(pow(10.,0.05*grid_gain),vmin=0,vmax=3.,extent=(0,360.,0.,90.))
    plt.scatter(gsm_var[:,1]*180./pi,90.-gsm_var[:,0]*180./pi,s=1,linewidth=0,c=gain_beam,vmin=0,vmax=3,cmap=cm.jet)
    plt.colorbar() 
    plt.xlim(0,360)
    plt.ylim(0,90)
#    pylab.xlabel('Azimuth (degrees)')
    plt.ylabel('Altitude (degrees)')
    plt.title('HIbiscus Beam (linear power)')
    plt.subplot(312)
#    pylab.imshow(grid_temp,vmin=0,vmax=2e4,extent=(0,360.,0.,90.))
    plt.scatter(gsm_var[:,1]*180./pi,90.-gsm_var[:,0]*180./pi,s=1,linewidth=0,c=gsm_array,vmin=0,vmax=2e4,cmap=cm.jet)
    plt.colorbar()
    plt.xlim(0,360)
    plt.ylim(0,90)
#    pylab.xlabel('Azimuth (degrees)') 
    plt.ylabel('Altitude (degrees)')
    plt.title('GSM Data (Kelvin)')
    plt.subplot(313) 
#    pylab.imshow(pow(10.,0.05*grid_gain)*grid_temp,vmin=0,vmax=2e4,extent=(0,360.,0.,90.))
    plt.scatter(gsm_var[:,1]*180./pi,90.-gsm_var[:,0]*180./pi,s=1,linewidth=0,c=full_beam,vmin=0,vmax=5e4,cmap=cm.jet)
    plt.colorbar() 
    plt.xlim(0,360)
    plt.ylim(0,90)
    plt.xlabel('Azimuth (degrees)')  
    plt.ylabel('Altitude (degrees)')
    plt.title('Expected Signal (Kelvin)')
    plt.subplots_adjust(hspace=0.4)
    plt.savefig(label,dpi=300)
    plt.clf()

#    grid_gain = []
#    grid_temp = []

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
