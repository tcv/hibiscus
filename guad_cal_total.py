import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf

#main_cal = 'Isla_Guadalupe_data_jun_2013/cal_data/'
#main_dir = $GUAD_DATA
#cal_dirs = os.listdir(main_dir)
cal_dir = '/home/tcv/lustre/cal_data/'
caldir = '/lustre/cal_data/'
cal_files = os.listdir(cal_dir)
#print cal_dirs

diff_time = []
Temp_gain = []
GainR = []
GainX = []
InR = []
InX = []
VnR = []
VnX = []
for i in range(1,15):
    date = str(i)
    if i<10:
        date = '0'+str(i)
    date_name = 'June'+date
    lim_Cal_time = loadtxt(cal_dir+date_name+'Cal_time.txt')
#    lim_Cal_time = array(lim_Cal_time)
    lim_Temp_gain=loadtxt(cal_dir+date_name+'TempGain_avg_take4.txt')
    lim_GainR = loadtxt(cal_dir+date_name+'Real_Gain_avg_take4.txt')
    lim_GainX=loadtxt(cal_dir+date_name+'Imag_Gain_avg_take4.txt')
    lim_InR=loadtxt(cal_dir+date_name+'Real_In_avg_take4.txt')
    lim_InX=loadtxt(cal_dir+date_name+'Imag_In_avg_take4.txt')
    lim_VnR=loadtxt(cal_dir+date_name+'Real_Vn_avg_take4.txt')
    lim_VnX=loadtxt(cal_dir+date_name+'Imag_Vn_avg_take4.txt')
    new_freq = loadtxt(cal_dir+date_name+'Cal_freq.txt')
    if len(lim_Temp_gain)<500:
	for i in range(0,len(lim_Cal_time)):
	    diff_time.append(lim_Cal_time[i])
            Temp_gain.append(lim_Temp_gain[i])
            GainR.append(lim_GainR[i])
            GainX.append(lim_GainX[i])
            InR.append(lim_InR[i])
            InX.append(lim_InX[i])
            VnR.append(lim_VnR[i])
            VnX.append(lim_VnX[i])
    else: 
#        for i in range(0,len(lim_Cal_time)):
        diff_time.append(lim_Cal_time)
        Temp_gain.append(lim_Temp_gain)
        GainR.append(lim_GainR) 
        GainX.append(lim_GainX)
        InR.append(lim_InR)
        InX.append(lim_InX)
        VnR.append(lim_VnR)
        VnX.append(lim_VnX)


diff_time = array(diff_time)
Temp_gain = array(Temp_gain)
GainR = array(GainR)
GainX = array(GainX)
InR = array(InR)
InX = array(InX)
VnR = array(VnR)
VnX = array(VnX)

single_GR = GainR[:,4000]
median_GR = ma.median(real(single_GR))
std_GR = ma.std(real(single_GR))
pylab.imshow(GainR,vmax=(median_GR+10*std_GR),vmin=0,aspect=100./len(GainR),extent=(40,140,len(GainR),0))
cbar = pylab.colorbar()
cbar.set_label('Real Gain')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Real Gain Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_Real_Gain_waterfall_take4',dpi=300)
pylab.clf()

single_GX = GainX[:,4000]
median_GX = ma.median(single_GX)
std_GX = ma.std(single_GX)
pylab.imshow(GainX,vmin=(median_GX-10*std_GX),vmax=0,aspect=100./len(GainX),extent=(40,140,len(GainX),0))
cbar = pylab.colorbar()
cbar.set_label('Imaginary Gain')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Imag Gain Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_Imag_Gain_waterfall_take4',dpi=300)
pylab.clf()

single_TG = Temp_gain[:,4000]
median_TG = ma.median(single_TG)
std_TG = ma.std(single_TG)
pylab.imshow(Temp_gain,vmax=(median_TG+100*std_TG),vmin=0,aspect=100./len(Temp_gain),extent=(40,140,len(Temp_gain),0))
cbar = pylab.colorbar()
cbar.set_label('Temperature Gain')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Temperature Gain Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_TempGain_waterfall_take4',dpi=300)
pylab.clf()

single_VR = VnR[:,4000]
median_VR = ma.median(single_VR)
std_VR = ma.std(single_VR)
pylab.imshow(VnR,vmax=(median_VR+10*std_VR),vmin=0,aspect=100./len(VnR),extent=(40,140,len(VnR),0))
cbar = pylab.colorbar()
cbar.set_label('Real Vn')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Real Vn Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_Real_Vn_waterfall_take4',dpi=300)
pylab.clf()

single_VX = VnX[:,4000]
median_VX = ma.median(single_VX)
std_VX = ma.std(single_VX)
pylab.imshow(VnX,vmax=(median_VX+10*std_VX),vmin=0,aspect=100./len(VnX),extent=(40,140,len(VnX),0))
cbar = pylab.colorbar()
cbar.set_label('Imag Vn')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Imag Vn Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_Imag_Vn_waterfall_take4',dpi=300)
pylab.clf()

single_IR = InR[:,4000]
median_IR = ma.median(single_IR)
std_IR = ma.std(single_IR)
pylab.imshow(InR,vmax=(median_IR+10*std_IR),vmin=0,aspect=100./len(InR),extent=(40,140,len(InR),0))
cbar = pylab.colorbar()
cbar.set_label('Real In')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Real In Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_Real_In_waterfall_take4',dpi=300)
pylab.clf()

single_IX = InX[:,4000]
median_IX = ma.median(single_IX)
std_IX = ma.std(single_IX)
pylab.imshow(InX,vmax=(median_IX+10*std_IX),vmin=0,aspect=100./len(InX),extent=(40,140,len(InX),0))
cbar = pylab.colorbar()
cbar.set_label('Imag In')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Cal Sample')
pylab.title('Imaginary In Data Variation over time')
pylab.savefig('/home/tcv/'+caldir+'Full_Imag_In_waterfall_take4',dpi=300)
pylab.clf()

print 'Number of calibration files:',len(Temp_gain)
single_TG = Temp_gain[:,4000]
median_TG = ma.median(single_TG)
print median_TG
good_ind = []
for i in range(0,len(single_TG)):
    if single_TG[i]<2*median_TG:
#    if single_TG[i]!=NaN:
         if single_TG[i]>0:
            good_ind.append(i)


print 'Good Files:',len(good_ind),'Total_Files:',len(single_TG)

#diff_time_sm = diff_time[good_ind]
#Temp_gain_sm = []
#GainR_sm = []
#GainX_sm = []
#InR_sm = []
#InX_sm = []
#VnR_sm = []
#VnX_sm = []
#funct = lambda p,x: p[0]*x+p[1]
#err = lambda p,x,y: funct(p,x)-y
#pG=[1e19,1e18]
#pE = [1e-7,1e-8]
#for i in range(0,len(Temp_gain[0])):
#    TG = itp.UnivariateSpline(diff_time[good_ind],Temp_gain[good_ind,i],s=0)
#    Temp_gain_sm.append(TG)
#    GR = itp.UnivariateSpline(diff_time[good_ind],GainR[good_ind,i],s=0)
#    GainR_sm.append(GR)
#    GX = itp.UnivariateSpline(diff_time[good_ind],GainX[good_ind,i],s=0)
#    GainX_sm.append(GX)
#    VR = itp.UnivariateSpline(diff_time[good_ind],VnR[good_ind,i],s=0)
#    VnR_sm.append(VR)
#    VX = itp.UnivariateSpline(diff_time[good_ind],VnX[good_ind,i],s=0)
#    VnX_sm.append(VX)
#    IR = itp.UnivariateSpline(diff_time[good_ind],InR[good_ind,i],s=0)
#    InR_sm.append(IR)
#    IX = itp.UnivariateSpline(diff_time[good_ind],InX[good_ind,i],s=0)
#    InX_sm.append(IX)

#print 'Shape of Smoothed Cal Data:',shape(InR_sm)
#print Temp_gain_sm[1000]
ant_time = arange(diff_time[0]-1,diff_time[-1]+1,0.01)
pylab.scatter(diff_time[good_ind],Temp_gain[good_ind,1000])
#pylab.scatter(ant_time,Temp_gain_sm[1000](ant_time))
#pylab.ylim(0,2e18)
pylab.savefig(cal_dir+'overall_cal_temp_take4',dpi=300)
pylab.clf()

mean_TG = ma.mean(Temp_gain[good_ind,:],axis=0)
mean_GR = ma.mean(GainR[good_ind,:],axis=0)
mean_GX = ma.mean(GainX[good_ind,:],axis=0)
mean_VR = ma.mean(VnR[good_ind,:],axis=0)
mean_VX = ma.mean(VnX[good_ind,:],axis=0)
mean_IR = ma.mean(InR[good_ind,:],axis=0)
mean_IX = ma.mean(InX[good_ind,:],axis=0)

print 'Number of Frequency Channels is:',shape(new_freq)
print 'Number of Mean Data Points is:',shape(mean_TG)
tot_mean = ma.mean(mean_TG)
std_mean = ma.std(mean_TG)
pylab.plot(new_freq,ones(len(mean_TG))*tot_mean)
pylab.plot(new_freq,ones(len(mean_TG))*(tot_mean+std_mean/2.))
pylab.plot(new_freq,ones(len(mean_TG))*(tot_mean-std_mean/2.))
pylab.plot(new_freq,mean_TG)
#pylab.ylim(0,5e12)
#pylab.ylim(0,5e11)
#pylab.ylim(0,2e18)
pylab.savefig(cal_dir+'overall_cal_temp_mean_take4',dpi=300)
pylab.clf()

savetxt(cal_dir+'Real_In_avg_take4.txt',mean_IR,delimiter = ' ')
savetxt(cal_dir+'Imag_In_avg_take4.txt',mean_IX,delimiter = ' ')
savetxt(cal_dir+'Real_Vn_avg_take4.txt',mean_VR,delimiter = ' ')
savetxt(cal_dir+'Imag_Vn_avg_take4.txt',mean_VX,delimiter = ' ')
savetxt(cal_dir+'Real_Gain_avg_take4.txt',mean_GR,delimiter =' ')
savetxt(cal_dir+'Imag_Gain_avg_take4.txt',mean_GX,delimiter =' ')
savetxt(cal_dir+'TempGain_avg_take4.txt',mean_TG,delimiter=' ')
savetxt(cal_dir+'Cal_time_take4.txt',diff_time,delimiter=' ')
savetxt(cal_dir+'Cal_freq_take4.txt',new_freq,delimiter=' ')
