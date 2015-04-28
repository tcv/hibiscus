from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import skrf as rf
import ephem as eph
import numpy.polynomial.polynomial as poly
import errno

def loadsingle(fname):
    """
    loads a single file information.
    Output is:
    file timestamp (in UTC Hours since the first of the month).
    file form (antenna, short, open, 50ohm)
    file data array
    voltage (from the header) - default is 0 Volts
    temperature (from the header) - default is 300 Kelvin
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
    Tamb=300.
    next_line = False
    time2 = 0.
    
    f = open(fname)
    for i in range(0,20):
        test = f.readline()
        if test.split(' ')[-1]=='V\n':
            volt = float(test.split(' ')[2])
        elif test.split(' ')[-1]=='celsius\n':
            Tamb = float(test.split(' ')[2])+273.15
        if next_line ==True:
            month = test.split('-')[1]
            end = test.split('-')[-1]
            day = end.split(' ')[0]
            tt = end.split(' ')[1]
            hr = tt.split(':')[0]
            minu = tt.split(':')[1]
            sec = tt.split(':')[2]
            time2 = float(day)*24.+float(hr)+float(minu)/60.+float(sec)/3600.

            next_line = False
        if len(test.split('  '))>1:

            if test.split('  ')[1]=='GMT Time:\n':
                next_line = True

    if abs(time-time2)>0.01:
        print 'Difference betweeen fname time and ftime: ',time-time2

    return time,form,file_data,mask_data,freq_data,volt,Tamb

def writesingle(fname, new_dir, new_data, append):
    """
    writes a single file in the same format 
    as the original file with updated data.
    takes as input the filename (with path) 
    for the original file, the new directory, 
    the new data and the addition to the file names 
    (if desired) - append='' gives no addition. 

    assumes that you want the same header as the original file. 
    """
    f = open(fname)
    file_name = fname.split('/')[-1]
    new_fname = new_dir+file_name.split('.')[0]+str(append)+'.dat'
    nf = open(new_fname,'w+')
    for line in f:
        if line[0]=='#':
            newline = str(line)
            nf.write(newline)
    
    for i in range(0,len(new_data)):
        dataline = str(new_data[i])+'\n'
        nf.write(dataline)
    
    nf.close()

    return

    

def flagging(data,freq,sigma_thres,linscale):
    """
    Flags data for RFI.
    Designed for a single time step scan.
    Uses a sigma threshold to flag out anything with
    RFI over a certain threshold.
   
    Also flags out NaNs, infs.
    Inputs are:
    data - linear input
    freq - can be any units
    sigma_thres - cutoff for bad data
    linscale - size of flattened window

    Output is flagging mask for input data array.
    """

    mask = zeros(len(data))
    nanmask = array(where(isnan(data))[0])
    mask[nanmask] = 1.0
    infmask = array(where(isinf(data))[0])
    mask[infmask] = 1.0
    scale = linscale
    for f in range(0, len(data)/scale-1):
        (Fa,Fb) = polyfit(freq[f*scale:(f+1)*scale],data[f*scale:(f+1)*scale],1)
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
        
    (Fa,Fb) = polyfit(freq[(f+1)*scale:-1],data[(f+1)*scale:-1],1)
    flat_data = data[(f+1)*scale:-1]/polyval([Fa,Fb],freq[(f+1)*scale:-1])
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


def threshold_flag(data,masked,freq,value):
    """
    Flags out RFI above a set cutoff.
    Cutoff has been defined using value as a multiplicative factor times the average polynomial guess (set to 1 at 70 MHz).
    """
    new_mask = zeros(len(data))
    for i in range(0,len(data)):
        if masked[i]==1.0:
            new_mask[i] = 1.0

    log_freq = log10(freq/70.)
    fit_params = [3.5186,-2.6205]
    dfit = 10**(poly.polyval(log_freq,fit_params))
    f70 = where(freq<=70.)[0][-1]
    subfit = dfit/dfit[f70]
    for i in range(0,len(data)):
        if new_mask[i]==0:
            data_fit = subfit[i]*data[f70]
            if data[i]>data_fit*value:
                new_mask[i] = 1.0

    return new_mask

def cal_flag(short_data,short_fit,masked,freq,cutoff):
    """
    Flags frequencies with bad short data compared to the fit.
    Note, this will flag frequencies outside the short fit band.
    """
    new_mask = zeros(len(masked))

    for i in range(0,len(masked)):
        if masked[i]==1.0:
            new_mask[i] = 1.0

    for i in range(0,len(masked)):
        if abs(short_data[i]-short_fit[i])>cutoff:
            new_mask[i] = 1.0
            
    return new_mask

def spike_flag(data,masked,freq,percent):
    """
    Flags out RFI spikes using a 11 bin filter
    Can be used with either time or freq
    percent is a percentage level cut (100 would be twice the 11 bin average)
    Needs to be applied to masked data.
    """
    new_mask = zeros(len(data))
    new_array = ma.array(data,mask=masked)
    new_comp = ma.compressed(new_array)
    freq_array = ma.array(freq,mask=masked)
    new_freq = ma.compressed(freq_array)
    for i in range(0,len(data)):
        if masked[i]==1.0:
            new_mask[i] = 1.0

    for i in range(5,len(new_comp)-5):
        group = new_comp[i-5]+new_comp[i-4]+new_comp[i-3]+new_comp[i-2]+new_comp[i-1]+new_comp[i]+new_comp[i+1]+new_comp[i+2]+new_comp[i+3]+new_comp[i+4]+new_comp[i+5]
        mean_group = group/11.
        if new_comp[i]/mean_group>=(1+percent/100.):
            comp_freq = new_freq[i]
            for j in range(0,len(freq)):
                if freq[j]==comp_freq:
                    index=j
            new_mask[index]= 1.0
        elif new_comp[i]/mean_group<=1/(1+percent/100.):
            comp_freq = new_freq[i]
            for j in range(0,len(freq)):
                if freq[j]==comp_freq:
                    index=j
            new_mask[index]= 1.0
   
    return new_mask

def timeflag(data,masked,time,sigma_thres,linscale):
    """
    Flags time data for RFI from single time steps 
    over a certain threshold
    Designed to be used for a single frequency dataset (axis is time).
    data - linear data to be flagged    
    masked - current mask for the data (skips already masked data)
    time - array of times used (analogous to freq array)
    sigma_thres - sigma threshold set for cutoff
    linscale - amount of time data to flatten

    Output is the updated mask
    """

    scale = linscale
    new_mask = zeros(len(data))
    for i in range(0,len(new_mask)):
        if masked[i]==1.0:
            new_mask[i]=1.0

    f=0
    if (len(data)/scale-1)<=1:
        (Fa,Fb) = polyfit(time,data,1)
        flat_data = data/polyval([Fa,Fb],time)
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
    Output is rebinned data with corresponding freq, mask arrays.
    """

    if binscale > 1:
        new_data = zeros(len(data)/binscale)
        new_mask = zeros(len(data)/binscale)
        new_freq = zeros(len(data)/binscale)
        f=0
        for f in range(0, len(new_data)-1):
            if len(masked[f*binscale:(f+1)*binscale])==sum(masked[f*binscale:(f+1)*binscale]):
                new_data[f] = 1.0
            else: 
                test_data = ma.array(data[f*binscale:(f+1)*binscale],mask=masked[f*binscale:(f+1)*binscale])
                test_data_con = ma.compressed(test_data)
                new_data[f] = ma.mean(test_data_con)
            if sum(masked[f*binscale:(f+1)*binscale])>=binscale/2.:
                new_mask[f] = 1.0
            new_freq[f] = ma.mean(freq[f*binscale:(f+1)*binscale])
        if len(masked[(f+1)*binscale:-1])==sum(masked[(f+1)*binscale:-1]):
            new_data[-1] = 1.0
        else:
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
        if sum(masked[:,f])==len(masked[:,f]):
            new_data[f] = 1.0
        else:
            masked_data = ma.array(data[:,f],mask=masked[:,f])
            compressed_data = ma.compressed(masked_data)
            new_data[f] = ma.mean(compressed_data)
        if sum(masked[:,f])>=len(data[0])/2.:
            new_mask[f] = 1.0
        

    return new_data,new_mask

def sidereal(time,idate,longitude,latitude):
    """
    Turns a time in UTC (defined since initial) into LST.
    time is in hours and idate is in 'YYYY/MM/DD'
    """
    initial = eph.date('2013/6/1')
    site = eph.Observer()
    site.lon = longitude
    site.lat = latitude
#    guad.lon = '-118.3'
#    guad.lat = '28.8833'
    
    single = eph.date(initial+time/24.)
    site.date = single
    single_time = site.sidereal_time()
    sidereal_hour = single_time*12./pi

    return sidereal_hour

def mkdir_p(path):
    """
    Create a directory with 'mkdir -p' functionality
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise

    return

