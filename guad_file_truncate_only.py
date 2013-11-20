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
outdir = '/lustre/tcv/truncated_data/'
#ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
#amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'
#cal_dir = '/home/tcv/lustre/cal_data_sept/'
#diff_time = loadtxt(cal_dir+'Cal_time_avg.txt')
#GainR = loadtxt(cal_dir+'Real_Gain_avg_fit.txt')
#GainX = loadtxt(cal_dir+'Imag_Gain_avg_fit.txt')
#InR = loadtxt(cal_dir+'Real_In_avg_fit.txt')
#InX = loadtxt(cal_dir+'Imag_In_avg_fit.txt')
#VnR = loadtxt(cal_dir+'Real_Vn_avg_fit.txt')
#VnX = loadtxt(cal_dir+'Imag_Vn_avg_fit.txt')
#cal_freq = loadtxt(cal_dir+'Cal_freq.txt')
#print shape(GainR)

#R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
#R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
#Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
#Eff_sm = fc.smooth(Effic,F_ant,0.01)
#R_ant_sm = fc.smooth(R_ant,F_ant,0)
#X_ant_sm = fc.smooth(X_ant,F_ant,0)
#R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
#X_amp_sm = fc.smooth(X_amp,F_amp,4e3)

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
#            print len(sub_data)
            width = 250.0/len(sub_data)
            freqs = arange(0,250.0,width)
            mask = zeros(len(sub_data))
            new_data,new_mask,new_freq = fc.truncate(sub_data,mask,freqs,40.,130.)
            new_data = array(new_data)
            new_data = 10**(new_data/10.)
            nandata = where(isnan(new_data))[0]
            for i in range(0,len(nandata)):
                new_data[nandata[i]] = 0.0
            infdata = where(isinf(new_data))[0]
            for i in range(0,len(infdata)):
                new_data[infdata[i]] = 0.0
            if volt>9.5:
                if volt<12.0:
                    fc.writesingle(filename,new_directory,new_data)
