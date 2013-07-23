from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf

def loadsingle(fname):
    """
    loads a single file information.
    Output is:
    file timestamp (in UTC Hours since the first of the month).
    file form (antenna, short, open, 50ohm)
    file data array
    voltage (from the header)
    if present may also have
    file mask array
    file frequency array
    """
    file_name = fname.split('/')[-1]
    end = file_name.split('-')[5]
    sec = end.split('_')[0]
    nofile = end.split('.')[0]
    form = nofile.split('_')[3]
    mask_data = []
    freq_data = []
    time =float(file_name.split('-')[2])*24.0+float(file_name.split('-')[3])+float(file_name.split('-')[4])/60+float(sec)/3600
    file_dump = loadtxt(fname)
    file_dump = array(file_dump)
    if len(file_dump)>5:
        file_data = file_dump
    elif len(file_dump)==2:
        file_data = file_dump[0]
        mask_data = file_dump[1]
    elif len(file_dump)==3:
        file_data = file_dump[0]
        mask_data = file_dump[1]
        freq_data = file_dump[2]
    else:
        file_data = [0.0,]
        mask_data = [0.0,]
        freq_data = [0.0,]

    volt = 0.0
    f = open(fname)
    test = f.readlines(1)
    for line in range(0,len(test)):
        if test[line].split('\t')[0]=='#  VIN:':
            lim = test[line].split('\t')[1]
            volt = float(lim.split('\n')[0])

    return time,form,file_data,mask_data,freq_data,volt

def imped(ant_data,cable_len):
    """
    Calculates impedence from the smith chart data.
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
    Can add a phase shift if needed.
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

def smooth(Eff_array,freq,s_factor):
    """
    Gives a smooth function version of a noisy array
    using an interpolation.
    Expects frequency values in MHz
    s_factor is:
    3e-20 for gain cal file
    0.3 for efficiency array
    0.475 for 50 ohm load efficiency array
    2e-19 for log period gain cal
    """
    func = itp.UnivariateSpline(freq,Eff_array,s=s_factor)
    
    return func

def effcal(Smooth_array, freq, sub_data):
    """
    Function for applying efficiency calibration function.
    Assumes that Efficiency only known above 25 MHz.
    Expects frequency in MHz.
    Expects data array to be linear.
    Output is also linear, sets all data below 25 MHz to zero.
    """
    corrected_data = zeros(len(sub_data))
    for f in range(0,len(sub_data)):
        if freq[f]>=25.0:
            corrected_data[f] = sub_data[f]/Smooth_array(freq[f])

    return corrected_data


def flagging(data,freq,sigma_thres,linscale):
    """
    Flags data for RFI.
    Designed for a single time step scan.
    Uses a sigma threshold to flag out anything with
    RFI over a certain threshold.
    Expects data to be linear for spline (s=1e-10). want to try something else.
    seems like using db data getting reasonable results for s = 1e4
    
    Also flags out NaNs, infs.

    Output is flagging mask for input data array.
    """
#    data = 10.**(data/10.)
    mask = zeros(len(data))
    nanmask = array(where(isnan(data))[0])
    mask[nanmask] = 1.0
    infmask = array(where(isinf(data))[0])
    mask[infmask] = 1.0
    scale = linscale
    for f in range(0, len(data)/scale-1):
 #       smooth = itp.UnivariateSpline(freq[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale])
	(Fa,Fb) = polyfit(freq[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale],1)
# 	smooth = itp.interp1d(freq[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale],'linear')
        flat_data = data[f*scale:(f+1)*scale]/polyval([Fa,Fb],freq[f*scale:(f+1)*scale])
        flat_sigma = ma.std(flat_data)
        flat_mean = ma.mean(flat_data)
        max_accept = 1.0+flat_sigma*sigma_thres
        min_accept = 1.0-flat_sigma*sigma_thres
        maxmask = array(where(flat_data>max_accept)[0])
        minmask = array(where(flat_data<min_accept)[0])
        maxmask = maxmask+f*scale
        minmask = minmask+f*scale
        mask[maxmask] = 1.0
        mask[minmask] = 1.0
        
#    smooth = itp.UnivariateSpline(freq[(f+1)*scale:-1],data[(f+1)*scale:-1])
#    smooth = itp.interp1d(freq[(f+1)*scale:-1],data[(f+1)*scale:-1],'linear')
    (Fa,Fb) = polyfit(freq[(f+1)*scale:-1],data[(f+1)*scale:-1],1)
    flat_data = data[(f+1)*scale:-1]/polyval([Fa,Fb],freq[(f+1)*scale:-1])
#    flat_data = data[(f+1)*scale:-1]/smooth(freq[(f+1)*scale:-1])
    flat_sigma = ma.std(flat_data)
    flat_mean = ma.mean(flat_data)
    max_accept = 1.0+flat_sigma*sigma_thres
    min_accept = 1.0-flat_sigma*sigma_thres
    maxmask = array(where(flat_data>max_accept)[0])
    minmask = array(where(flat_data<min_accept)[0])
    maxmask = maxmask+(f+1)*scale
    minmask = minmask+(f+1)*scale
    mask[maxmask] = 1.0
    mask[minmask] = 1.0
    
    return mask

def spike_flag(data,percent):
    """
    Flags out RFI spikes using a 11 bin filter
    Can be used with either time or freq
    percent is a percentage level cut (100 would be twice the 11 bin average)
    """
    mask = zeros(len(data))
    for i in range(5,len(data)-5):
        group = data[i-5]+data[i-4]+data[i-3]+data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2]+data[i+3]+data[i+4]+data[i+5]
        mean_group = group/11.
        if data[i]/mean_group>(1+percent/100):
            mask[i]= 1.0
        elif data[i]/mean_group<1/(1+percent/100):
            mask[i]=1.0
    
    return mask

def timeflag(data,masked,time,sigma_thres,linscale):
    """
    Flags time data for RFI
    Designed to be used for a single frequency dataset (axis is time).
    Uses a sigma threshold to flag out outliers.
    linscale sets the amount of time data to scale over (~20)
    """
#    data = 10.**(data/10.)
    scale = linscale
    new_mask = zeros(len(data))
    for i in range(0,len(new_mask)):
        if masked[i]==1.0:
            new_mask[i]=1.0

    f=0
    if (len(data)/scale-1)<=1:
#	smooth = itp.interp1d(time,data,'linear')
#        smooth = itp.UnivariateSpline(time,data)
        (Fa,Fb) = polyfit(time,data,1)
	flat_data = data/polyval([Fa,Fb],time)
#        flat_data = data/smooth(time)
        flat_sigma = ma.std(flat_data)
        flat_mean = ma.mean(flat_data)
        max_accept = 1.0+flat_sigma*sigma_thres
        min_accept = 1.0-flat_sigma*sigma_thres
        maxmask = array(where(flat_data>max_accept)[0])
        minmask = array(where(flat_data<min_accept)[0])
        maxmask = maxmask
        minmask = minmask
        new_mask[maxmask] = 1.0
        new_mask[minmask] = 1.0
    elif(len(data)/scale-1)>1:
        for f in range(0, len(data)/scale-1):
#	    smooth = itp.interp1d(time[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale],'linear')
#            smooth = itp.UnivariateSpline(time[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale])
#            flat_data = data[f*scale:(f+1)*scale]/smooth(time[f*scale:(f+1)*scale])
	    (Fa,Fb) = polyfit(time[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale],1)
	    flat_data = data[f*scale:(f+1)*scale]/polyval([Fa,Fb],time[f*scale:(f+1)*scale])
            flat_sigma = ma.std(flat_data)
            flat_mean = ma.mean(flat_data)
            max_accept = 1.0+flat_sigma*sigma_thres
            min_accept = 1.0-flat_sigma*sigma_thres
            maxmask = array(where(flat_data>max_accept)[0])
            minmask = array(where(flat_data<min_accept)[0])
            maxmask = maxmask+f*scale
            minmask = minmask+f*scale
            new_mask[maxmask] = 1.0
            new_mask[minmask] = 1.0
        
#	smooth = itp.interp1d(time[(f+1)*scale:-1],data[(f+1)*scale:-1])
#        smooth = itp.UnivariateSpline(time[(f+1)*scale:-1],data[(f+1)*scale:-1])
#        flat_data = data[(f+1)*scale:-1]/smooth(time[(f+1)*scale:-1])
	(Fa,Fb) = polyfit(time[(f+1)*scale:-1],data[(f+1)*scale:-1],1)
	flat_data = data[(f+1)*scale:-1]/polyval([Fa,Fb],time[(f+1)*scale:-1])
        flat_sigma = ma.std(flat_data)
        flat_mean = ma.mean(flat_data)
        max_accept = 1.0+flat_sigma*sigma_thres
        min_accept = 1.0-flat_sigma*sigma_thres
        maxmask = array(where(flat_data>max_accept)[0])
        minmask = array(where(flat_data<min_accept)[0])
        maxmask = maxmask+(f+1)*scale
        minmask = minmask+(f+1)*scale
        new_mask[maxmask] = 1.0
        new_mask[minmask] = 1.0

    return new_mask

    
def rebin(data,masked,freq,binscale):
    """
    Rebins data to coarser frequency resolution.
    Assumes that the input is the data after flagging, mask,
    corresponding freq array and a binscale (number of bins to merge).
    Uses median rather than mean to avoid single freq outliers
    Output is rebinned data with corresponding freq, mask arrays.
    """
    if binscale > 1:
        new_data = zeros(len(data)/binscale)
        new_mask = zeros(len(data)/binscale)
        new_freq = zeros(len(data)/binscale)
        f=0
        for f in range(0, len(new_data)-1):
            test_data = ma.array(data[f*binscale:(f+1)*binscale],mask=masked[f*binscale:(f+1)*binscale])
            test_data_con = ma.compressed(test_data)
            new_data[f] = ma.mean(test_data_con)
            if sum(masked[f*binscale:(f+1)*binscale])>=binscale/2.:
                new_mask[f] = 1.0
            new_freq[f] = ma.mean(freq[f*binscale:(f+1)*binscale])
        test_data = ma.array(data[(f+1)*binscale:-1],mask=masked[(f+1)*binscale:-1])
        test_data_con = ma.compressed(test_data) 
        new_data[-1] = ma.mean(test_data_con)
        if sum(masked[(f+1)*binscale:-1])>=1.:
            new_mask[-1] = 1.0
        new_freq[-1] = ma.mean(freq[(f+1)*binscale:-1])
    elif binscale == 1:
        new_data = data
        new_mask = masked
        new_freq = freq
    
    return new_data,new_mask,new_freq

def truncate(data,mask,freq,minfreq,maxfreq):
    """
    Removes the data that is outside a given range to shrink array size.
    Expects freq, minfreq, maxfreq to be in same units. 
    """
    new_data = []
    new_freq = []
    new_mask = []
    for i in range(0,len(freq)):
        if freq[i]>minfreq:
            if freq[i]<maxfreq:
                new_data.append(data[i])
                new_freq.append(freq[i])
                new_mask.append(mask[i])
    new_data = array(new_data)
    new_mask = array(new_mask)
    new_freq = array(new_freq)

    return new_data,new_mask,new_freq

def timerebin(data,masked):
    """
    Rebins chunk of time data to single dataset
    Assumes that the input is a two dimensional array with corresponding mask
    Output is single dataset with corresponding mask
    """
    new_data = zeros(len(data[0]))
    new_mask = zeros(len(data[0]))
    data = array(data)
    masked = array(masked)
    
    for f in range(0,len(data[0])):
        masked_data = ma.array(data[:,f],mask=masked[:,f])
        compressed_data = ma.compressed(masked_data)
        new_data[f] = ma.mean(compressed_data)
        if sum(masked[:,f])>=len(data[0])/2.:
            new_mask[f] = 1.0
        

    return new_data,new_mask

def waterfallplot(data, Tscale, time, freq, name,filename):
    """
    Creates a waterfall plot of data.
    Expects data to be in time,freq array.
    """
    diff = 0
    setting = time[0]
    for i in range(1,len(time)):
        if time[i]>setting:
            diff+=time[i]-setting
            setting = time[i]
        elif time[i]<setting:
            diff+=(24.0-setting)+time[i]
            setting = time[i]

    print diff        
    asp = (freq[-1]-freq[0])/diff

    pylab.imshow(data,vmin=Tscale[0],vmax=Tscale[1],aspect=asp,
                 extent=(freq[0],freq[-1],diff,0.0))
    pylab.colorbar()
    pylab.title(name)
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Duration (Hours), Start=%0.2f UTC' %(time[0]%24.0))
    pylab.savefig(filename,dpi=200)
    pylab.clf()

    return


def noise_calc(load,short,term,noise,Zamp,freq):
    """
    Calculates the Vn and In for a given 50 Ohm, 100 Ohm, Noise and Short dataset.
    Assumes 50 Ohm, 100 Ohm, Noise and short are in dB, freq is in MHz
    """
    Z100 = 100.*e**(2*pi*freq*400*1e-6*1j)
    Z50 = 50.*ones(len(freq))
    Vsh = sqrt(2*Z50*10**(short/10.))
    V50 = sqrt(2*Z50*10**(load/10.))
    V100 = sqrt(2*Z50*10**(load/10.))
#    Vsh = sqrt(absolute(Zamp)*10**(short/10.))
#    V50 = sqrt(absolute(Zamp)*10**(load/10.))
#    V100 = sqrt(absolute(Zamp)*10**(term/10.))
    VJF = sqrt(4*1.381e-23*300*Z50*(freq[1]-freq[0])*1e6)*ones(len(freq))
#    VJO = sqrt(4*1.381e-23*300*absolute(Z100)*(freq[1]-freq[0])*1e6)*ones(len(freq))
    VJO = sqrt(2.)*VJF
    Gain = (Z100*(V50-Zamp*Vsh/(Zamp+Z50))*(Zamp+Z50)-Z50*(V100-Zamp*Vsh/(Zamp+Z100))*(Zamp+Z100))/(Zamp*VJF*(Z100-sqrt(2.)*Z50))
    In = (V50-Zamp*Vsh/(Zamp+Z50))*(Zamp+Z50)/(Gain*Zamp*Z50)-VJF/Z50 
##    In = VJF*((Zamp+Z50)*sqrt(2.)*(V50-Vsh)/((Zamp+Z100)*(V100-Vsh))-1)/(Z50-Z100*(Zamp+Z50)*(V50-Vsh)/((Zamp+Z100)*(V100-Vsh)))
##    Gain = (V50-Vsh)*(Zamp+Z50)/(Zamp*(VJF+In*Z50))
#    Gain = ((Zamp+Z50)*Z100*V50+Zamp*(Z50-Z100)*Vsh-(Zamp+Z100)*V100*Z50)/(Zamp*(VJF*Z100-VJO*Z50))
    Vn = Vsh/Gain
#    In = ((Zamp+Z50)*V50/(Zamp*Z50)-(Vsh+Gain*VJF)/Z50)/Gain
#    Vnoise = sqrt(absolute(Zamp)*10**(noise/10.))
    Vnoise = sqrt(2*Z50*10**(noise/10.))
##    Vns = Vnoise/Gain - Zamp*(VJF+Vn)/(Zamp+Z50)-In*Zamp*Z50/(Zamp+Z50)
    Vm50 = Zamp*(VJF+Vn)/(Zamp+Z50)+In*Zamp*Z50/(Zamp+Z50)
    Vns = sqrt(((Zamp+Z50)*Vnoise/Zamp/Gain)**2-Vm50**2)
#    Pns = absolute(Vns**2/Zamp)
    Pnoise = abs(Vnoise)**2
    P50 = abs(Vm50*Gain)**2
    Psh = abs(Vn*Gain)**2
#    P50 = absolute(V50**2/(Zamp*Gain**2))
#    Vm50 = Zamp*(VJF+Vn)/(Zamp+Z50)+In*Zamp*Z50/(Zamp+Z50)
#    Ve50 = V50/Gain
#    Vns = sqrt((Vnoise/Gain)**2-(Vm50)**2)
#    Tns = 300*(abs(Vnoise/Gain)**2-abs(Vm50)**2)/abs(Vm50)**2   
#    Tns = 300*(abs(Vnoise)**2)*(Zamp+Z50)/Zamp/Gain/abs(Vm50)**2
    Pns = abs(Vns*Gain)**2/(2*Z50)
    Tns = 300.*Pnoise/(P50-Psh)
    gtemp = Tns/Pns

    return Vn,In,Gain,gtemp

def noise_corr(data,Vn,In,freq,Zamp,Zant,Gain,Temp):
    """
    Subtracts the Vn and In from a dataset and converts to Vsky.
    Assumes data is in dB and freq is in MHz
    """
    Pant = 10**(data/10.)
    Vn1 = Vn
    In1 = In
#    Vant = sqrt(absolute(Zamp)*Pant)
    Vant = sqrt(2*50.*Pant)
    Vsky = Vant*(Zamp+Zant)/Zamp/Gain-Vn1-In1*Zant
    Psky = absolute(Vsky*Gain)**2*real(Zant)/(2*absolute(Zant)**2)
#    Psky = absolute(Vsky**2/Zant)
    Tsky = Temp*Psky

    return Tsky

