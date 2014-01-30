import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import sys

maindir = '/lustre/anat/Guadalupe_data/'
outdir = '/lustre/tcv/calibrated_data_test/'
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'
cal_dir = '/home/tcv/lustre/cal_data_sept/'
#diff_time = loadtxt(cal_dir+'Cal_time_avg.txt')
GainR = loadtxt(cal_dir+'Real_Gain_avg_fit.txt')
GainX = loadtxt(cal_dir+'Imag_Gain_avg_fit.txt')
InR = loadtxt(cal_dir+'Real_In_avg_fit.txt')
InX = loadtxt(cal_dir+'Imag_In_avg_fit.txt')
VnR = loadtxt(cal_dir+'Real_Vn_avg_fit.txt')
VnX = loadtxt(cal_dir+'Imag_Vn_avg_fit.txt')
cal_freq = loadtxt(cal_dir+'Cal_freq.txt')
print shape(GainR)

R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)
R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)

directories = os.listdir(maindir)

date_ind = sys.argv[1]

if int(date_ind)<15:
    direct = 'June'+date_ind+'_day_to_night'
    directory = maindir+direct+'/'
    new_directory = outdir+direct+'/'
    dirlist = os.listdir(directory)
    for fname in dirlist:
        if len(fname.split('_'))>=3:
            filename = directory+fname
            time,form,sub_data,mask,freq,volt,temp = fc.loadsingle(filename)
            print len(sub_data)
            width = 250.0/len(sub_data)
            freqs = arange(0,250.0,width)
            Z_ant = R_ant_sm(freqs)+1j*X_ant_sm(freqs)
            Z_amp = R_amp_sm(freqs)+1j*X_amp_sm(freqs)
            Z50 = 50.*ones(len(freqs))
            Z100 = 100.*exp(2*pi*freqs*400*1e-6*1j)
            cal_ind = 0
            Tsky = zeros(len(sub_data))

#            for k in range(0,len(diff_time)):
#                if (time-diff_time[k])>0.1:
#                    cal_ind = k
            TGain_s = 0.0 
#            Vn_s = VnR[cal_ind]+1j*VnX[cal_ind]
#            In_s = InR[cal_ind]+1j*InX[cal_ind]
#            Gain_s = GainR[cal_ind]+1j*GainX[cal_ind]
            Vn_s = []
            In_s = []
            Gain_s = []
            for i in range(0,len(freqs)):
                if freqs[i]<cal_freq[0]:
                    Vn_s.append(0.0)
                    In_s.append(0.0)
                    Gain_s.append(1.0)
                elif freqs[i]>cal_freq[-1]:
                    Vn_s.append(0.0)
                    In_s.append(0.0)
                    Gain_s.append(1.0)
                else:
                    for j in range(0,len(cal_freq)):
                        if freqs[i]==cal_freq[j]:
                            Vn_s.append(VnR[j]+1j*VnX[j])
                            In_s.append(InR[j]+1j*InX[j])
                            Gain_s.append(GainR[j]+1j*GainX[j])
            if form=='antenna':
                for j in range(0,len(sub_data)):
                    TGain_s = 0.0
                    Tsky[j] = fc.noise_corr(sub_data[j],Vn_s[j],In_s[j],width,Z_amp[j],Z_ant[j],Gain_s[j],TGain_s)
                corr_data = fc.effcal(Eff_sm,freqs,Tsky)
            elif form=='50ohm':
                lin_data = 10.**(sub_data/10.)
                Gt = (1/(4*1.381e-23*abs(Z_amp)*width*1e6))
                Tsky = 2*50.*lin_data*Gt/absolute(Gain_s)**2
            elif form=='open':
                lin_data = 10.**(sub_data/10.)
                Gt = 1/(4*1.381e-23*abs(Z_amp)*width*1e6)
                Tsky = 2*50.*lin_data*Gt/absolute(Gain_s)**2
            elif form=='short':
                lin_data = 10.**(sub_data/10.)
                Gt = 1/(4*1.381e-23*abs(Z_amp)*width*1e6)
                Tsky = 2*50.*lin_data*Gt/absolute(Gain_s)**2
            elif form=='noise':
                lin_data = 10.**(sub_data/10.)
                Gt = 1/(4*1.381e-23*abs(Z_amp)*width*1e6)
                Tsky = 2*50.*lin_data*Gt/absolute(Gain_s)**2
            fc.writesingle(filename,new_directory,Tsky)
