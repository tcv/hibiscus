import matplotlib
matplotlib.use('Agg')
from numpy import *
import matplotlib.pyplot as pylab
#from pylab import *
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
from scipy import linalg
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import ephem as eph

#data_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/'
data_dir = '/home/tcv/lustre/processed_data_take4/'
result_dir = '/home/tcv/lustre/data_plots_take4/'
gsm_raw_data = loadtxt('/home/tcv/guad_extras/gsm_guadalupe.dat')

data = os.listdir(data_dir)

processed_data = []
processed_time = []
processed_freq = []

dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
date_files = []

gsm_freq = arange(50,111,1)
gsm_time = arange(0,24,24./288)
gsm_data = zeros((len(gsm_time),len(gsm_freq)))
for t in range(0,len(gsm_time)):
    for f in range(0,len(gsm_freq)):
        gsm_data[t,f] = gsm_raw_data[t*61+f,2]

for i in range(0,len(dates)):
    curr_date = 'June'+dates[i]
    print 'Date being processed is:',curr_date
    single_data = loadtxt(data_dir+curr_date+'_processed_data_avgcal.txt')
    single_data = array(single_data)
    for k in range(0,len(single_data)):
	processed_data.append(single_data[k])
    single_time = loadtxt(data_dir+curr_date+'_processed_time_avgcal.txt')
    for j in range(0,len(single_time)):
        processed_time.append(single_time[j])
    if len(processed_freq)<1:
        processed_freq = loadtxt(data_dir+curr_date+'_processed_freq_avgcal.txt')
    print 'Current Size of Processed Data is:',shape(processed_data)    
	    
#    for file in data:
#	curr_date = 'June'+dates[i]
#	if file.split('-')[0]==curr_date:
#	    if file.split('.')[-1]=='txt':
#	        date_files.append(file)
#    print 'Number of Files in that day is:',len(date_files)/3
#    num_files = len(date_files)/3
#    for j in range(0,num_files):
#	full_date = curr_date+'-'+str(j)
#	if full_date != 'June09-5':
#      	    single_data = loadtxt(data_dir+full_date+'_processed_data_avgcal.txt')
# 	    single_data = array(single_data)
#	    if len(single_data)<550:
#	        for k in range(0,len(single_data)):
#		    processed_data.append(single_data[k])
#	    else:
#	        processed_data.append(single_data)
#	    single_time = loadtxt(data_dir+full_date+'_processed_time_avgcal.txt')
#	    processed_time.append(single_time)
#	    if len(processed_freq)<1:
#	        processed_freq = loadtxt(data_dir+full_date+'_processed_freq_avgcal.txt')
#	    print 'Current status of processed data:',shape(processed_data)
#    date_files = []


#print 'Full Shape of Processed Data is:',shape(processed_data)
processed_data = array(processed_data)
processed_time = array(processed_time)
print 'Full Shape of Processed Data is:',shape(processed_data)
print 'Shape of Processed Time is:',shape(processed_time)
print where(isnan(processed_data))[0],where(isinf(processed_data))[0]

# Test SVD:
U,S,Vh = linalg.svd(processed_data[:,101:184])
#print U.shape,S.shape,Vh.shape
Sm = linalg.diagsvd(S,len(processed_data),len(processed_data[0,101:184]))
S0 = S
S0[0] = 0. 
Sm0 = linalg.diagsvd(S0,len(processed_data),len(processed_data[0,101:184]))
P1 = dot(U,dot(Sm0,Vh))
S0[1] = 0.
Sm1 = linalg.diagsvd(S0,len(processed_data),len(processed_data[0,101:184]))
P2 =  dot(U,dot(Sm1,Vh))
S0[2] = 0.
Sm2 = linalg.diagsvd(S0,len(processed_data),len(processed_data[0,101:184]))
P3 =  dot(U,dot(Sm2,Vh))

prelim_mask = zeros((len(processed_data),len(processed_data[0])))

for i in range(0,len(prelim_mask)):
    spike_mask = fc.spike_flag(processed_data[i],100)
    for j in range(0,len(prelim_mask[0])):
	if spike_mask[j] == 1.0:
	    prelim_mask[i,j] = 1.0
	    processed_data[i,j] = 0.0
        if processed_data[i,j]==0.0:
            prelim_mask[i,j] = 1.0

mean_data = []
for i in range(0,len(processed_data[0])):
    single_freq = ma.array(processed_data[:,i],mask=prelim_mask[:,i])
    single_compress = ma.compressed(single_freq)
    single_mean = ma.mean(single_compress)
    mean_data.append(single_mean)
nandata = where(isnan(mean_data))
nandata = array(nandata)    
for i in range(0,len(nandata[0])):
    index=nandata[0,i]
    mean_data[index]=0.0

fitfunc = lambda p,x: p[0]*x**(-2.5)+p[1]
errfunc = lambda p,x,y: fitfunc(p,x)-y
p0 = [10.,10.]
p,success = opt.leastsq(errfunc,p0,args=(processed_freq,mean_data))
p1,success = opt.leastsq(errfunc,p,args=(processed_freq[81:205],mean_data[81:205]))
p2,success = opt.leastsq(errfunc,p1,args=(processed_freq[101:184],mean_data[101:184]))

initial = eph.date('2013/6/1')
ephem_dates = []
for i in range(0,len(processed_time)):
    ephem_dates.append(eph.date(initial+processed_time[i]/24.))
    

guad = eph.Observer()
guad.lon = '-118.3'
guad.lat = '28.8833'

ephem_sidereal = []
for i in range(0,len(processed_time)):
    guad.date = ephem_dates[i]
    test = guad.sidereal_time()
    ephem_sidereal.append(test)

sidereal_hour = []
for i in range(0,len(processed_time)):
    sidereal_hour.append(ephem_sidereal[i]*12./pi)

sidereal_sort = argsort(sidereal_hour)
sidereal_sort_time = sort(sidereal_hour)
sidereal_data = zeros((len(processed_data),len(processed_data[0])))
for i in range(0,len(sidereal_sort)):
    sidereal_data[i] = processed_data[sidereal_sort[i]]

stack_time = gsm_time
stack_data = zeros((len(stack_time),len(processed_freq)))
for i in range(0,len(stack_time)):
    sub_data = []
    for j in range(0,len(sidereal_sort_time)):
        if abs(stack_time[i]-sidereal_sort_time[j])<1./24.:
            sub_data.append(sidereal_data[j])
    stack_data[i] = ma.mean(sub_data,axis=0)
nandata = where(isnan(stack_data))
nandata = array(nandata)
for i in range(0,len(nandata[0])):
    index=nandata[0,i]
    stack_data[index]=zeros(len(stack_data[0]))
infdata = where(isinf(stack_data))
infdata = array(infdata)
for i in range(0,len(infdata[0])):
   index=infdata[0,i]
   stack_data[index]=zeros(len(stack_data[0]))


Us,Ss,Vhs = linalg.svd(stack_data[:,101:184])
Ssm = linalg.diagsvd(Ss,len(stack_data),len(stack_data[0,101:184]))
Ss0 = Ss
Ss0[0] = 0.
Ssm0 = linalg.diagsvd(Ss0,len(stack_data),len(stack_data[0,101:184]))
Ps1 = dot(Us,dot(Ssm0,Vhs))
Ss0[1] = 0.
Ssm1 = linalg.diagsvd(Ss0,len(stack_data),len(stack_data[0,101:184]))
Ps2 = dot(Us,dot(Ssm1,Vhs))
Ss0[2] = 0.
Ssm2 = linalg.diagsvd(Ss0,len(stack_data),len(stack_data[0,101:184]))
Ps3 = dot(Us,dot(Ssm2,Vhs))

full_time = arange(processed_time[0],processed_time[-1],0.05)
full_data = zeros((len(full_time),len(processed_freq)))
for i in range(0,len(full_data)):
    for j in range(0,len(processed_time)):
        if abs(full_time[i]-processed_time[j])<=0.1:
            full_data[i]=processed_data[j]

Uf,Sf,Vhf = linalg.svd(full_data[:,101:184])
Sfm = linalg.diagsvd(Sf,len(full_data),len(full_data[0,101:184]))
Sf0 = Sf
Sf0[0] = 0.
Sfm0 = linalg.diagsvd(Sf0,len(full_data),len(full_data[0,101:184]))
Pf1 = dot(Uf,dot(Sfm0,Vhf))
Sf0[1] = 0.
Sfm1 = linalg.diagsvd(Sf0,len(full_data),len(full_data[0,101:184]))
Pf2 = dot(Uf,dot(Sfm1,Vhf))
Sf0[2] = 0.
Sfm2 = linalg.diagsvd(Sf0,len(full_data),len(full_data[0,101:184]))
Pf3 = dot(Uf,dot(Sfm2,Vhf))

full_time_sidereal = []
for t in range(0,len(full_time)):
    single_date = eph.date(initial+full_time[t]/24.)
    guad.date = single_date
    single_time = guad.sidereal_time()
    single_hour = single_time*12./pi
    full_time_sidereal.append(single_hour)

#pylab.plot(full_time_sidereal)
#pylab.savefig('test_time')
#pylab.clf()

#print gsm_freq[0], processed_freq[0]
print where(processed_freq<65)[0][-1]
print where(processed_freq>85)[0][0]

full_gsm_data = zeros((len(full_time),len(processed_freq)))
for t in range(0,len(full_time)):
    min_time = 0
    max_time = 0
    for tg in range(0,len(gsm_time)):
	if (float(full_time_sidereal[t])-float(gsm_time[tg]))>0:
            min_time = tg 
            max_time = tg+1
        elif float(full_time_sidereal[t])==float(gsm_time[tg]):
            min_time = tg
            max_time = tg        
    if float(full_time_sidereal[t])>float(gsm_time[-1]):
        min_time = -1
        max_time = -1
    for f in range(0,len(processed_freq)):
	min_freq = 0
        max_freq = 0
        for fg in range(0,len(gsm_freq)):
            if (float(processed_freq[f])-float(gsm_freq[fg]))>0:
                min_freq = fg
                max_freq = fg+1
            elif float(processed_freq[f])==float(gsm_freq[fg]):
                min_freq = fg
                max_freq = fg
        if float(processed_freq[f])<=float(gsm_freq[0]):
            min_freq = 0
            max_freq = 0
        elif float(processed_freq[f])>=float(gsm_freq[-1]):
            min_freq = -1
            max_freq = -1
	if full_data[t,f]>0.0:
#            if min_time == max_time:
#                if min_freq == max_freq:
#                    gsm_diff[t,f] = full_data[t,f]-gsm_data[min_time,min_freq]
#                else:
#                    (Ta,Tb) = polyfit(gsm_freq[min_freq:max_freq],gsm_data[min_time,min_freq:max_freq],1)
#                    gsm_diff[t,f] = full_data[t,f]-polyval([Ta,Tb],processed_freq[f])
#            else:
#                if min_freq == max_freq:
#                    (Ta,Tb) = polyfit(gsm_time[min_time:max_time],gsm_data[min_time:max_time,min_freq],1)
#                    gsm_diff[t,f] = full_data[t,f]-polyval([Ta,Tb],full_time[t]) 
#                else:
#		    (Tamin,Tbmin) = polyfit(gsm_freq[min_freq:max_freq],gsm_data[min_time,min_freq:max_freq],1)
#                    (Tamax,Tbmax) = polyfit(gsm_freq[min_freq:max_freq],gsm_data[max_time,min_freq:max_freq],1)
#                    tdata = [polyval([Tamin,Tbmin],processed_freq[f]),polyval([Tamax,Tbmax],processed_freq[f])]
#                    (Ta,Tb) = polyfit(gsm_time[min_time:max_time],tdata,1)
#                    gsm_diff[t,f] = full_data[t,f]-polyval([Ta,Tb],full_time[t])
            if processed_freq[f]<50.:
                full_gsm_data[t,f]=0.0
            elif processed_freq[f]>110.:
                full_gsm_data[t,f]=0.0 
            else:
                full_gsm_data[t,f] = gsm_data[min_time,min_freq]
#    print min_time     

pylab.imshow(full_data,vmin=0,vmax=2e4,aspect=100./(full_time[-1]-full_time[0]),extent=(40,140,full_time[-1],full_time[0]))
cbar = pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_avgcal',dpi=300)
pylab.clf()

pylab.imshow(full_gsm_data,vmin=0,vmax=2e4,aspect=100./(full_time[-1]-full_time[0]),extent=(40,140,full_time[-1],full_time[0]))
cbar = pylab.colorbar() 
cbar.set_label('Temperature (Kelvin)') 
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
pylab.title('GSM Model Data for full dataset time/freq')
#pylab.title('Full Dataset Waterfall Plot - GSM Model Data') 
pylab.savefig(result_dir+'full_waterfall_gsm_avgcal',dpi=300) 
pylab.clf() 

pylab.imshow(full_data[:,40:287]-full_gsm_data[:,40:287],vmin=-2e3,vmax=2e3,aspect=60./(full_time[-1]-full_time[0]),extent=(50.,110.,full_time[-1],full_time[0]))
cbar = pylab.colorbar()
cbar.set_label('Temperature Difference (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset - GSM Model Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_data-gsm_avgcal',dpi=300)
pylab.clf()

pylab.scatter(processed_time,processed_data[:,120],c='b',edgecolor='b',s=3)
pylab.xlabel('Time (Hours since Midnight June 01 GMT)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Single Frequency (%0.1f MHz) Time Variation' %processed_freq[120])
pylab.grid()
pylab.xlim(processed_time[0],processed_time[-1])
pylab.ylim(0,5e3)
pylab.savefig(result_dir+'single_freq_time_var_avgcal',dpi=300)
pylab.clf()

pylab.scatter(sidereal_hour,processed_data[:,120],c='b',edgecolor='b',s=3)
pylab.xlabel('Time (Sidereal Hour)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Single Frequency (%0.1f MHz) Time Variation Folded by Day' %processed_freq[120])
pylab.grid()
pylab.xlim(0,24)
pylab.ylim(0,5e3)
pylab.savefig(result_dir+'single_freq_time_var_fold_avgcal',dpi=300)
pylab.clf()

pylab.scatter(processed_freq,mean_data,c='b',edgecolor='b',s=3)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Mean Data Spectrum')
pylab.grid()
pylab.xlim(40,140)
pylab.ylim(0,2e4)
pylab.savefig(result_dir+'mean_freq_spectrum_avgcal',dpi=300)
pylab.plot(processed_freq,fitfunc(p,processed_freq),c='g')
pylab.plot(processed_freq,fitfunc(p1,processed_freq),c='r')
pylab.plot(processed_freq,fitfunc(p2,processed_freq),c='c')
pylab.legend(('Mean Data','Full Data Fit','Fit to 60-90 MHz','Fit to 65-85 MHz'))
pylab.title('Mean Data Spectrum with -2.5 Power Law Fits')
pylab.savefig(result_dir+'mean_freq_spectrum_avgcal_fit',dpi=300) 
pylab.clf()

pylab.scatter(sidereal_hour,processed_data[:,100],c='b',edgecolor='b',s=3)
pylab.scatter(sidereal_hour,processed_data[:,125],c='g',edgecolor='g',s=3)
pylab.scatter(sidereal_hour,processed_data[:,150],c='r',edgecolor='r',s=3)
pylab.scatter(sidereal_hour,processed_data[:,175],c='c',edgecolor='c',s=3)
pylab.scatter(sidereal_hour,processed_data[:,200],c='m',edgecolor='m',s=3)
pylab.xlabel('Time (Sidereal Hour)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Multiple Frequency Time Variation Folded by Day')
pylab.grid()
pylab.xlim(0,24)
pylab.ylim(0,1e4)
pylab.legend(('%0.1f MHz' %processed_freq[100],'%0.1f MHz' %processed_freq[125],'%0.1f MHz' %processed_freq[150],'%0.1f MHz' %processed_freq[175],'%0.1f MHz' %processed_freq[200]))
pylab.savefig(result_dir+'multi_freq_time_var_fold_avgcal',dpi=300)
pylab.clf()


### SVD Plots
pylab.imshow(P1,vmin=-1e2,vmax=3e2,aspect=20./(processed_time[-1]-processed_time[0]),extent=(65,85,processed_time[-1],processed_time[0]))
cbar= pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
pylab.title('Full Dataset Waterfall Plot')
pylab.savefig(result_dir+'Unfill_waterfall_1SVDrm_avgcal_65-85MHz',dpi=300) 
pylab.clf() 

pylab.imshow(P2,vmin=-1e2,vmax=3e2,aspect=20./(processed_time[-1]-processed_time[0]),extent=(65,85,processed_time[-1],processed_time[0]))
cbar= pylab.colorbar() 
cbar.set_label('Temperature (Kelvin)') 
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')  
pylab.title('Full Dataset Waterfall Plot') 
pylab.savefig(result_dir+'Unfill_waterfall_2SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()  

pylab.imshow(P3,vmin=-1e2,vmax=3e2,aspect=20./(processed_time[-1]-processed_time[0]),extent=(65,85,processed_time[-1],processed_time[0]))
cbar= pylab.colorbar()  
cbar.set_label('Temperature (Kelvin)')  
pylab.xlabel('Frequency (MHz)')   
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot')  
pylab.savefig(result_dir+'Unfill_waterfall_3SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()   

pylab.imshow(processed_data[:,101:184]-dot(U,dot(Sm,Vh)),vmin=-1,vmax=1,aspect=20./(processed_time[-1]-processed_time[0]),extent=(65,85,processed_time[-1],processed_time[0]))
cbar= pylab.colorbar()  
cbar.set_label('Temperature (Kelvin)')  
pylab.xlabel('Frequency (MHz)')   
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot')  
pylab.savefig(result_dir+'Unfill_waterfall_data-AllSVD_avgcal_65-85MHz',dpi=300)
pylab.clf()   

pylab.imshow(Pf1,vmin=-1e2,vmax=3e2,aspect=20./(full_time[-1]-full_time[0]),extent=(65,85,full_time[-1],full_time[0]))
cbar= pylab.colorbar() 
cbar.set_label('Temperature (Kelvin)') 
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot') 
pylab.savefig(result_dir+'full_waterfall_1SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()  
 
pylab.imshow(Pf2,vmin=-1e2,vmax=3e2,aspect=20./(full_time[-1]-full_time[0]),extent=(65,85,full_time[-1],full_time[0]))
cbar= pylab.colorbar()  
cbar.set_label('Temperature (Kelvin)')  
pylab.xlabel('Frequency (MHz)')   
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot')  
pylab.savefig(result_dir+'full_waterfall_2SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()   

pylab.imshow(Pf3,vmin=-1e2,vmax=3e2,aspect=20./(full_time[-1]-full_time[0]),extent=(65,85,full_time[-1],full_time[0]))
cbar= pylab.colorbar()   
cbar.set_label('Temperature (Kelvin)')   
pylab.xlabel('Frequency (MHz)')    
pylab.ylabel('Time (Hours since Midnight June 01 GMT)') 
pylab.title('Full Dataset Waterfall Plot')   
pylab.savefig(result_dir+'full_waterfall_3SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()    
 
pylab.imshow(full_data[:,101:184]-dot(Uf,dot(Sfm,Vhf)),vmin=-1,vmax=1,aspect=20./(full_time[-1]-full_time[0]),extent=(65,85,full_time[-1],full_time[0]))
cbar= pylab.colorbar()   
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')  
pylab.ylabel('Time (Hours since Midnight June 01 GMT)')
pylab.title('Full Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_data-AllSVD_avgcal_65-85MHz',dpi=300)
pylab.clf()

pylab.imshow(Ps1,vmin=-5e2,vmax=5e2,aspect=20./24.,extent=(65,85,24,0))
cbar= pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Sidereal Hour)')
pylab.title('Stack Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_stack_1SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()
  
pylab.imshow(Ps2,vmin=-50,vmax=50,aspect=20./24.,extent=(65,85,24,0))
cbar= pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Sidereal Hour)')
pylab.title('Stack Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_stack_2SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()
 
pylab.imshow(Ps3,vmin=-50,vmax=50,aspect=20./24.,extent=(65,85,24,0))
cbar= pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Sidereal Hour)')
pylab.title('Stack Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_stack_3SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()
  
pylab.imshow(stack_data[:,101:184]-dot(Us,dot(Ssm,Vhs)),vmin=-1,vmax=1,aspect=20./24.,extent=(65,85,24,0))
cbar= pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Sidereal Hour)')
pylab.title('Stack Dataset Waterfall Plot')
pylab.savefig(result_dir+'full_waterfall_stack_data-AllSVD_avgcal_65-85MHz',dpi=300)
pylab.clf()

svd1_mean_data = []
svd2_mean_data = []
svd3_mean_data = []
for i in range(0,len(P1[0])):
#    single_freq = ma.array(processed_data[:,i],mask=prelim_mask[:,i])
    single_freq_1svd = ma.array(P1[:,i],mask=prelim_mask[:,i])
    single_freq_2svd = ma.array(P2[:,i],mask=prelim_mask[:,i]) 
    single_freq_3svd = ma.array(P3[:,i],mask=prelim_mask[:,i])  
    single_compress_1svd = ma.compressed(single_freq_1svd)
    single_compress_2svd = ma.compressed(single_freq_2svd)
    single_compress_3svd = ma.compressed(single_freq_3svd)
#    single_compress = ma.compressed(single_freq)
#    single_mean = ma.mean(single_compress)
    single_mean_1svd = ma.mean(single_compress_1svd)
    svd1_mean_data.append(single_mean_1svd)
    single_mean_2svd = ma.mean(single_compress_2svd)
    svd2_mean_data.append(single_mean_2svd)
    single_mean_3svd = ma.mean(single_compress_3svd)
    svd3_mean_data.append(single_mean_3svd)
nandata1 = where(isnan(svd1_mean_data))
nandata1 = array(nandata1)
for i in range(0,len(nandata1[0])):
    index=nandata1[0,i]
    svd1_mean_data[index]=0.0
nandata2 = where(isnan(svd2_mean_data))
nandata2 = array(nandata2)
for i in range(0,len(nandata2[0])):
    index=nandata2[0,i]
    svd2_mean_data[index]=0.0
nandata3 = where(isnan(svd3_mean_data))
nandata3 = array(nandata3)
for i in range(0,len(nandata3[0])):
    index=nandata3[0,i]
    svd3_mean_data[index]=0.0

pylab.scatter(processed_freq[101:184],svd1_mean_data,c='b',edgecolor='b',s=3)
pylab.scatter(processed_freq[101:184],svd2_mean_data,c='g',edgecolor='g',s=3)
pylab.scatter(processed_freq[101:184],svd3_mean_data,c='r',edgecolor='r',s=3)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Mean Data Spectrum')
pylab.grid()
pylab.xlim(65,85)
pylab.ylim(-50,50)
pylab.legend(('1 SVD Mode Rm','2 SVD Mode Rm','3 SVD Mode Rm'))
pylab.savefig(result_dir+'mean_freq_spectrum_SVDrm_avgcal_65-85MHz',dpi=300)
pylab.clf()

pylab.imshow(stack_data,vmin=0,vmax=2e4,aspect=100./24.,extent=(40.,140.,24.,0.))
cbar = pylab.colorbar()
cbar.set_label('Temperature (Kelvin)')
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Time (Sidereal Hour)')
pylab.title('Full Dataset Waterfall Folded on Sidereal Hour')
pylab.savefig(result_dir+'full_waterfall_folded_avgcal',dpi=300)
pylab.clf()

svd1_mean_data = []
svd2_mean_data = []
svd3_mean_data = []
for i in range(0,len(Ps1[0])):
    single_mean_1svd = ma.mean(Ps1[:,i])
    svd1_mean_data.append(single_mean_1svd)
    single_mean_2svd = ma.mean(Ps2[:,i])
    svd2_mean_data.append(single_mean_2svd)
    single_mean_3svd = ma.mean(Ps3[:,i])
    svd3_mean_data.append(single_mean_3svd)
nandata1 = where(isnan(svd1_mean_data))
nandata1 = array(nandata1)
for i in range(0,len(nandata1[0])):
    index=nandata1[0,i]
    svd1_mean_data[index]=0.0
nandata2 = where(isnan(svd2_mean_data))
nandata2 = array(nandata2)
for i in range(0,len(nandata2[0])):
    index=nandata2[0,i]
    svd2_mean_data[index]=0.0
nandata3 = where(isnan(svd3_mean_data))
nandata3 = array(nandata3)
for i in range(0,len(nandata3[0])):
    index=nandata3[0,i]
    svd3_mean_data[index]=0.0

pylab.scatter(processed_freq[101:184],svd1_mean_data,c='b',edgecolor='b',s=3)
pylab.scatter(processed_freq[101:184],svd2_mean_data,c='g',edgecolor='g',s=3)
pylab.scatter(processed_freq[101:184],svd3_mean_data,c='r',edgecolor='r',s=3)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Mean Data Spectrum')
pylab.grid()
pylab.xlim(65,85)
pylab.ylim(-50,50)
pylab.legend(('1 SVD Mode Rm','2 SVD Mode Rm','3 SVD Mode Rm'))
pylab.savefig(result_dir+'mean_freq_spectrum_stack_SVDrm_avgcal_65-85MHz',dpi=300)
pylab.ylim(-5,5)
pylab.savefig(result_dir+'mean_freq_spectrum_stack_SVDrm_avgcal_65-85MHz_zoom',dpi=300)
pylab.clf()


