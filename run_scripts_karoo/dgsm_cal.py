"""
Module to do dGSM calibration.  
Assumes that the data has already been put into truncated npy arrays.
Also assumes that flagging masks have already been calculated, but will run without masks as well.

Takes 4 inputs (input directory, output directory, earliest time file day and hour in that directory)

Note that it also has inputs for the supplemental data needed to run the code. 
"""
import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
import os
import skrf as rf
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('../../hibiscus'))
import file_funcs as ff
import gsm_funcs as gf
import cal_funcs as cf
import eff_funcs as ef
import ephem as eph
import errno
import numpy.polynomial.polynomial as poly

#Setting input and output directories 
indir = sys.argv[1]
outdir = sys.argv[2]
supdir= '../../supplemental_data/'

#Data positions and times
gbt_idate = '2013/5/1'
gbt_lat = '38.433'
gbt_lon = '-79.8397'

guad_idate = '2013/6/1'
guad_lat = '28.8833'
guad_lon = '-118.3'

karoo_idate = '2015/4/1'
karoo_lat = '-30.7216'
karoo_lon = '21.4109'

#S parameter files
ant_s11_file =supdir+'SCIHI-070-v3.S1P'
amp_s_file = supdir+'WEA101_AMP_2013-04-04.s2p'
ant_s11_file_old = supdir+'ANT70MHZ_CMU.s1p'

# Efficiency calculation
# 'old' antenna data taken at CMU, other data taken in the Karoo
R_ant,X_ant,F_ant = ef.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = ef.imped_skrf(amp_s_file,0.0)
R_ant_old,X_ant_old,F_ant_old = ef.imped_skrf(ant_s11_file_old,0.0)
Z_ant = sqrt(R_ant**2+X_ant**2)
Z_ant_sm = itp.UnivariateSpline(F_ant,Z_ant,s=0.01)
Effic = ef.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = itp.UnivariateSpline(F_ant,Effic,s=0.01)
Effic_old = ef.effic(R_ant_old,X_ant_old,F_ant_old,R_amp,X_amp,F_amp)
Eff_sm_old = itp.UnivariateSpline(F_ant_old,Effic_old,s=0.01)

#Best fit parameters to short data collected in the Karoo.
fitshort = [9.43299312e-9,-1.16197388e-10,4.31005321e-13]

print sys.argv[1]
ant = sys.argv[1].split('_')[-1]
gsm_freqs = arange(50,111,1)
#Load file for gsm data with 70 MHz antenna.
if ant=='100/':
    gsm_data = numpy.load(supdir+'gsm_data_100_Karoo_test/gsm_data_full_100_Karoo.npy')
    gsm_times = numpy.load(supdir+'gsm_data_100_Karoo_test/gsm_sid_time_full_100_Karoo.npy')
    f70_gsm = where(gsm_freqs<=85.)[0][-1]
elif ant=='70/':
    gsm_data = numpy.load(supdir+'gsm_data_70_Karoo_test/gsm_data_full_70_Karoo.npy')
    gsm_times = numpy.load(supdir+'gsm_data_70_Karoo_test/gsm_sid_time_full_70_Karoo.npy')
    f70_gsm = where(gsm_freqs<=70.)[0][-1]

max_gsm = where(gsm_data[:,f70_gsm]==max(gsm_data[:,f70_gsm]))[0]
print max_gsm
print len(gsm_times)

directories = os.listdir(indir)

#Load files for full day of data. 
#Set to max length, then truncate. Max length is 1 dataset every 3 seconds for 24 hours. 
#Freq length set by previous experience, and is 11796 (truncated, but not rebinned).
data = zeros((28800,11796))
mask = zeros((28800,11796))
times = zeros(28800)
day_idate = int(sys.argv[3])*24+int(sys.argv[4])
dint = 0
mint = 0
tint = 0
fscale = 32
tscale = 32

for hour in range(0,24):
    for direct in directories:
         curr_day =int(direct.split('-')[2])
         end = direct.split('-')[-1]
         cur_hour = int(end.split('_')[0])
         cur_date = curr_day*24+cur_hour
         if (cur_date == day_idate+hour):
             if direct.split('_')[-1]=='antenna.npy':
                 single = numpy.load(indir+direct) 
                 if len(single)<2000:
                     freqs = arange(40.,130.,90./len(single[0]))
                     for time in range(0,len(single)):
                         data[dint] = single[time]
                         dint +=1
                 print 'Hour being processed is April ',curr_day,':',cur_hour

             elif direct.split('_')[-2] =='ant':
                 single = numpy.load(indir+direct)
                 for time in range(0,len(single)):
                     times[tint] = single[time]
                     tint+=1

             elif direct.split('_')[-1]=='mask.npy':
                 single = numpy.load(indir+direct)
                 for time in range(0,len(single)):
                     mask[mint] = single[time]
                     mint+=1


print 'Percent of Data Flagged from Frequency Masking: ',100.*sum(mask)/(len(mask)*len(mask[0]))
data = data[0:dint]
mask = mask[0:mint]
times = times[0:tint]

#Short fit for the given frequency range.
short_data = poly.polyval(freqs,fitshort)
st,sf,short_data,sm,sf,sv,ste = ff.loadsingle(supdir+'2015-04-05-00-06-26_Ch_2_noisetrunc.dat')

#Convert data times to sidereal times and sort.
sid_times = zeros(len(times))
for i in range(0,len(times)):
    single_sid = ff.sidereal(times[i],karoo_idate,karoo_lon,karoo_lat)
    sid_times[i] = single_sid


sortind = argsort(sid_times)
sorttime = zeros(len(sid_times))
sortdata = zeros((len(sid_times),len(data[0])))
sortmask = zeros((len(sid_times),len(data[0])))
for i in range(0,len(sid_times)):
    sorttime[i] = sid_times[sortind[i]]
    sortdata[i] = data[sortind[i]]
    sortmask[i] = mask[sortind[i]]
#    sortdata[i] = (data[sortind[i]]-short_data)/Eff_sm(freqs)
#    sortdata[i] = data[sortind[i]]-short_data


#Adding in Time Masking/threshold masking
for i in range(0,len(freqs)):
    new_mask = ff.timeflag(sortdata[:,i],sortmask[:,i],sorttime,3.,tscale)
    sortmask[:,i] = new_mask

#for i in range(0,len(sortdata)):
#    new_mask = ff.threshold_flag(sortdata[i],sortmask[i],freqs,75.)
#    sortmask[i] = new_mask

percent_masked = 100.*sum(sortmask)/(len(sortmask)*len(sortmask[0]))
print 'Percentage of Masked Data from Frequency and Time Masking: ',percent_masked

#Compress data to same time intervals as gsm data.
stack_data,stack_mask = cf.match_binning(gsm_times,freqs,sorttime,sortdata,sortmask)

#Correct for time error in data using time of max signal.
#70 MHz used as proxy for full dataset. 
#Re-sort based upon the adjustment. 

if ant=='100/':
    f70_gsm = where(gsm_freqs<=85.)[0][-1]
    f70_data = where(freqs<=85.)[0][-1]
elif ant=='70/':
    f70_gsm = where(gsm_freqs<=70.)[0][-1]
    f70_data = where(freqs<=70.)[0][-1]


data_70 = zeros(len(times))
for i in range(0,len(times)):
    data_70[i] = data[i][f70_data]

print stack_data[:,f70_data]
max_stack = where(stack_data[:,f70_data]==max(stack_data[:,f70_data]))[0]
print max_stack
#max_gsm = where(gsm_data[:,f70_gsm]==max(gsm_data[:,f70_gsm]))[0]
time_diff = gsm_times[max_gsm]-gsm_times[max_stack]
ind_diff = max_gsm-max_stack
adj_times = zeros(len(gsm_times))
for i in range(0,len(gsm_times)):
#    print gsm_times[i]
    adj_times[i] = (gsm_times[i]+time_diff)%24.

sortstack_ind = argsort(adj_times)
sortstack = zeros((len(adj_times),len(stack_data[0])))
sortstackmask=zeros((len(adj_times),len(stack_data[0])))
for i in range(0,len(sortstack_ind)):
    sortstack[i] = stack_data[sortstack_ind[i]]
    sortstackmask[i] = stack_mask[sortstack_ind[i]]

#Limit the gsm data to times where we have antenna data.
#Also extend the gsm data to the frequencies where we have data. 
lim_stack,lim_mask,lim_gsm,lim_time = cf.lim_bin(freqs,sortstack,sortstackmask,gsm_freqs,gsm_data,gsm_times)

#Calculate the time mean for the day.
mean_data,mean_mask = cf.time_mean(lim_stack,lim_mask)
gsm_mask = zeros((len(lim_gsm),len(lim_gsm[0])))
mean_gsm, mean_gmask = cf.time_mean(lim_gsm,gsm_mask)

#Generate mean subtracted arrays for calibration 
lim_ms_stack= zeros((len(lim_stack),len(stack_data[0])))
for i in range(0,len(lim_stack)):
    lim_ms_stack[i] = lim_stack[i]-mean_data
lim_ms_gsm = zeros((len(lim_stack),len(freqs)))
for i in range(0,len(lim_stack)):
    lim_ms_gsm[i] = lim_gsm[i]-mean_gsm

lim_ms_stack = array(lim_ms_stack)
lim_ms_gsm = array(lim_ms_gsm)

#Use least squares to calculate the calibration factor at each frequency
Kdgsm = zeros(len(freqs))
K0 = [6.e10,]
for i in range(0,len(freqs)):
    if len(lim_mask[:,i])==sum(lim_mask[:,i]):
        Kdgsm[i] = 1.e10
        print 'Fitting failed at frequency: ', freqs[i]
    else:
        Ki = cf.gain_calc(lim_ms_stack[:,i],lim_mask[:,i],lim_ms_gsm[:,i],K0[0])
        K0[0] = Ki
        Kdgsm[i] = Ki

#Plot the comparison of gsm to calibrated data at 70 MHz over time. 
pylab.plot(gsm_times,gsm_data[:,f70_gsm]-mean_gsm[f70_data]*ones(len(gsm_data)),label='gsm',c='k')
pylab.ylim(-3e3,3e3)
pylab.xlim(0,24)
pylab.xlabel('Sidereal Time (Hours)')
pylab.ylabel('Temperature (Kelvin)')
pylab.grid()
pylab.scatter(lim_time,Kdgsm[f70_data]*lim_ms_stack[:,f70_data],label='cal data',c='g',edgecolor='g',s=3)
pylab.savefig(outdir+'gsm_cal_test_'+sys.argv[3]+'_'+str(int(freqs[f70_data]))+'_MHz.png',dpi=300)
pylab.clf()


#Plot the time mean data using both the calculate calibration factor and a flat calibration factor (used as the initial input to the calibration). 
mean_array = ma.array(mean_data,mask=mean_mask)
mean_comp = ma.compressed(mean_array)
freq_comp = ma.compressed(ma.array(freqs,mask=mean_mask))
K_comp = ma.compressed(ma.array(Kdgsm,mask=mean_mask))
pylab.scatter(freq_comp,mean_comp*K_comp,c='b',edgecolor='b',s=3)
pylab.scatter(freq_comp,mean_comp*6e10,c='g',edgecolor='g',s=3)
pylab.ylim(0,1e4)
pylab.xlim(40,130)
pylab.grid()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.savefig(outdir+'gsm_cal_test_'+sys.argv[3]+'_mean.png',dpi=300)
pylab.clf()

if ant=='70/':
   fit_min = 50.
   fit_max = 90.
elif ant=='100/':
   fit_min = 80.
   fit_max = 110.
Kfit,Kparams = cf.poly_fore(mean_data*Kdgsm,mean_mask,freqs,fit_min,fit_max,2,ones(len(mean_data)))
print Kparams
print ma.min(mean_data*Kdgsm-Kfit)
print ma.max(mean_data*Kdgsm-Kfit)
print ma.mean(mean_data*Kdgsm-Kfit)
print ma.std(mean_data*Kdgsm-Kfit)

pylab.scatter(freqs,mean_data*Kdgsm-Kfit,c='b',edgecolor='b',s=3)
pylab.xlim(40,130)
pylab.ylim(-100,100)
pylab.grid()
pylab.xlabel('Frequency (MHz)')
pylab.savefig(outdir+'gsm_cal_test_'+sys.argv[3]+'_mean_fit.png',dpi=300)
pylab.clf()

cal_data = zeros((len(sortstack),len(sortstack[0])))
cal_data_masked = zeros((len(sortstack),len(sortstack[0])))
for t in range(0,len(cal_data)):
    cal_data[t] = sortstack[t]*Kdgsm
    cal_data_masked[t] = cal_data[t]
    for f in range(0,len(freqs)):
        if sortstackmask[t,f]==1:
            cal_data_masked[t,f] = 0.

pylab.imshow(cal_data,vmin=0,vmax=1e4,aspect=90./24.,extent=(freqs[0],freqs[-1],gsm_times[-1],gsm_times[0]))
pylab.colorbar()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.title('GSM Calibrated Temperature (Kelvin)')
pylab.savefig(outdir+'gsm_cal_test_'+sys.argv[3]+'_waterfall.png',dpi=300)
pylab.clf()

pylab.imshow(cal_data_masked,vmin=0,vmax=1e4,aspect=90./24.,extent=(freqs[0],freqs[-1],gsm_times[-1],gsm_times[0]))
pylab.colorbar()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.title('Masked GSM Cal Temperature (Kelvin)')
pylab.savefig(outdir+'gsm_cal_test_'+sys.argv[3]+'_masked_waterfall.png',dpi=300)
pylab.clf()

#numpy.save(outdir+'gsm_cal_data_Apr_'+sys.argv[3]+'_70MHz_ant.npy',cal_data_masked)
#numpy.save(outdir+'gsm_cal_times_Apr_'+sys.argv[3]+'_70MHz_ant.npy',gsm_times)
#numpy.save(outdir+'gsm_cal_values_Apr_'+sys.argv[3]+'_70MHz_ant.npy',Kdgsm)
if ant=='100/':
    numpy.save(outdir+'gsm_cal_data_Apr_'+sys.argv[3]+'_100MHz_ant.npy',cal_data_masked)
    numpy.save(outdir+'gsm_cal_times_Apr_'+sys.argv[3]+'_100MHz_ant.npy',gsm_times)
    numpy.save(outdir+'gsm_cal_values_Apr_'+sys.argv[3]+'_100MHz_ant.npy',Kdgsm)
elif ant=='70/':
    numpy.save(outdir+'gsm_cal_data_Apr_'+sys.argv[3]+'_70MHz_ant.npy',cal_data_masked)
    numpy.save(outdir+'gsm_cal_times_Apr_'+sys.argv[3]+'_70MHz_ant.npy',gsm_times)
    numpy.save(outdir+'gsm_cal_values_Apr_'+sys.argv[3]+'_70MHz_ant.npy',Kdgsm)


