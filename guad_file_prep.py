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

#maindir = 'Isla_Guadalupe_data_jun_2013/June15_day_to_night/'
#day_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/June15/'
maindir = '/lustre/anat/Guadalupe_data/'
daydir = '/lustre/data_arrays/'
rebin_ant = []
ant_time = []
ant_mask = []
ant_volt = []
binscale = 1
timescale = 1
ant_limit = []
mask_limit = []
time_limit = []
dataset = 0
directories = os.listdir(maindir)
new_freq = [0,1]

date_ind = sys.argv[1]

#dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
#for date_ind in range(9,10):
#for direct in directories:
if int(date_ind)<15:
    direct = 'June'+date_ind+'_day_to_night'
#    if direct.split('-')[0]=='2013':
    if direct.split('_')[-1]=='night':
        directory = maindir+direct+'/'
        print directory
        date = direct.split('_')[0]
#        print shape(rebin_ant)
#        if len(direct.split('_'))<2:
	dirlist = os.listdir(directory)
        for fname in dirlist:
	    if len(fname.split('_'))>=3:
	        filename=directory+fname
                time,form,sub_data,mask,freq,volt = fc.loadsingle(filename)
                width = 250.0/len(sub_data)
                freqs = arange(0,250.0,width)
                if len(sub_data)>1:
                    if form=='antenna':
#                mask = zeros(len(sub_data))
                        masked_sub = 10**(sub_data/10.)
                        mask = fc.flagging(masked_sub,freqs,3.,30)
                        excess_mask = fc.spike_flag(masked_sub,100)
                        for i in range(0,len(excess_mask)):
                            if excess_mask[i]==1.0:
                                mask[i] = 1.0
                        new_data,new_mask,new_freq = fc.rebin(sub_data,mask,freqs,binscale)
                        new_data,new_mask,new_freq = fc.truncate(new_data,new_mask,new_freq,40.,140.)
                        rebin_ant.append(new_data)
                        ant_time.append(time)
                        ant_mask.append(new_mask)
                        ant_volt.append(float(volt))
#                        print volt
	time_sort = argsort(ant_time)
#	print time_sort
	rebin_ant_sort = []
	ant_time_sort = [] 
	ant_mask_sort = []
        ant_volt_sort = []
	for t in range(0,len(rebin_ant)):
	    rebin_ant_sort.append(rebin_ant[time_sort[t]])
	    ant_time_sort.append(ant_time[time_sort[t]])
	    ant_mask_sort.append(ant_mask[time_sort[t]])
            ant_volt_sort.append(ant_volt[time_sort[t]])
	
	i=-1
	print 'Number of datasets in Day:',len(rebin_ant)
        print shape(ant_volt_sort)
	for i in range(0,len(rebin_ant)/500):
	    set = str(i)
            savetxt('/home/tcv/'+daydir+date+'-'+set+'_antenna.txt',rebin_ant_sort[i*500:(i+1)*500],delimiter = ' ')
            savetxt('/home/tcv/'+daydir+date+'-'+set+'_time.txt',ant_time_sort[i*500:(i+1)*500],delimiter=' ')
            savetxt('/home/tcv/'+daydir+date+'-'+set+'_mask.txt',ant_mask_sort[i*500:(i+1)*500],delimiter=' ')
            savetxt('/home/tcv/'+daydir+date+'-'+set+'_freq.txt',new_freq,delimiter=' ')
            savetxt('/home/tcv/'+daydir+date+'-'+set+'_volt.txt',ant_volt_sort[i*500:(i+1)*500],delimiter=' ')
	set = str(i+1)
        savetxt('/home/tcv/'+daydir+date+'-'+set+'_antenna.txt',rebin_ant_sort[(i+1)*500:-1],delimiter = ' ')
        savetxt('/home/tcv/'+daydir+date+'-'+set+'_time.txt',ant_time_sort[(i+1)*500:-1],delimiter=' ')
        savetxt('/home/tcv/'+daydir+date+'-'+set+'_mask.txt',ant_mask_sort[(i+1)*500:-1],delimiter=' ')
        savetxt('/home/tcv/'+daydir+date+'-'+set+'_freq.txt',new_freq,delimiter=' ') 
        savetxt('/home/tcv/'+daydir+date+'-'+set+'_volt.txt',ant_volt_sort[(i+1)*500:-1],delimiter=' ')

	print 'Number of Files in Day:',int(set)+1
        rebin_ant = []
        ant_time = []
        ant_mask = []
        ant_volt = []
#                        if len(time_limit)<1:
#                            time_test = time
#                        else:
#                            time_test = time_limit[-1]
#                        if time-time_test>0.1:
#                            timebin_data,timebin_mask = fc.timerebin(ant_limit,mask_limit)
#                            rebin_ant.append(timebin_data)
#                            ant_time.append(ma.mean(time_limit[0:-2]))
#                            ant_mask.append(timebin_mask)
#                            ant_limit = []
#                            mask_limit = []
#                            time_limit = []
#                            dataset = 0                            
#                        time_limit.append(time)
#                        ant_limit.append(new_data)
#                        mask_limit.append(new_mask)
#                        dataset = dataset + 1
#                        if dataset == timescale:
#                            timebin_data,timebin_mask = fc.timerebin(ant_limit,mask_limit)
#                            rebin_ant.append(timebin_data)
#                            ant_time.append(ma.mean(time_limit))
#                            ant_mask.append(timebin_mask)
#                            ant_limit = []
#                            mask_limit = []
#                            time_limit = []
#                            dataset = 0
