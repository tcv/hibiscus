import matplotlib
matplotlib.use('Agg')
from numpy import *
import pylab
#from pylab import *
#import matplotlib.pyplot as plt
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
#from scipy import optimize
import os
import data_analysis_funcs as fc
import skrf as rf
import sys
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#cal_dir = 'Isla_Guadalupe_data_jun_2013/cal_data/'
cal_dir = '/home/tcv/lustre/cal_data/'
#data_dir = 'Isla_Guadalupe_data_jun_2013/data_arrays/June15/'
data_dir = '/home/tcv/lustre/data_arrays/'
out_dir = '/home/tcv/lustre/processed_data_noeff/'
#out_dir = '/home/tcv/lustre/processed_data/'
#ant_s11_file = 'Isla_Guadalupe_data_jun_2013/ANT_3_average.s1p'
ant_s11_file = '/home/tcv/guad_extras/ANT_3_average.s1p'
#amp_s_file = 'Isla_Guadalupe_data_jun_2013/WEA101_AMP_2013-04-04.s2p'
amp_s_file = '/home/tcv/guad_extras/WEA101_AMP_2013-04-04.s2p'

diff_time = loadtxt(cal_dir+'Cal_time_take4.txt')
Temp_gain = loadtxt(cal_dir+'TempGain_fit_take4.txt')
GainR = loadtxt(cal_dir+'Real_Gain_fit_take4.txt')
GainX = loadtxt(cal_dir+'Imag_Gain_fit_take4.txt')
InR = loadtxt(cal_dir+'Real_In_fit_take4.txt')
InX = loadtxt(cal_dir+'Imag_In_fit_take4.txt')
VnR = loadtxt(cal_dir+'Real_Vn_fit_take4.txt')
VnX = loadtxt(cal_dir+'Imag_Vn_fit_take4.txt')
new_freq = loadtxt(cal_dir+'Cal_freq_take4.txt')

#single_TG = Temp_gain[:,1000]
#mean_TG = ma.median(single_TG)
#good_ind = []
#for i in range(0,len(single_TG)):
#    if single_TG[i]<2*mean_TG:
#    if single_TG[i]!=NaN:
#         if single_TG[i]>0:
#            good_ind.append(i)

##good_ind = [0,]
#print len(good_ind),len(single_TG)
#print good_ind

##Temp_gain_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
##GainR_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
##GainX_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
##InR_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
##InX_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
##VnR_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
##VnX_sm = zeros((2,len(good_ind),len(Temp_gain[0])))
#diff_time_sm = diff_time[good_ind]
##diff_time_sm = diff_time
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
##    for j in range(0,len(good_ind)-1):
##        xind = [diff_time[good_ind[j]],diff_time[good_ind[j+1]]]
##        yG = [Temp_gain[good_ind[j],i],Temp_gain[good_ind[j+1],i]]
##        TG,success = opt.leastsq(err,pG[:],args=(xind,yG))
##        pG = TG
#    TG = itp.UnivariateSpline(diff_time[good_ind],Temp_gain[good_ind,i],s=0)
#    Temp_gain_sm.append(TG)
##        Temp_gain_sm[:,j,i] = TG
##        yGR = [GainR[good_ind[j],i],GainR[good_ind[j+1],i]]
##        GR,success = opt.leastsq(err,pE[:],args=(xind,yGR))
##        pE = GR
#    GR = itp.UnivariateSpline(diff_time[good_ind],GainR[good_ind,i],s=0)
#    GainR_sm.append(GR)
##        GainR_sm[:,j,i] = GR
##        yGX =[GainX[good_ind[j],i],GainX[good_ind[j+1],i]]
##        GX,success = opt.leastsq(err,pE[:],args(xind,yGX))
#    GX = itp.UnivariateSpline(diff_time[good_ind],GainX[good_ind,i],s=0)
#    GainX_sm.append(GX)
##        GainX_sm[:,j,i] = GX
##        yVR = [VnR[good_ind[j],i],VnR[good_ind[j+1],i]]
##        VR,success = opt.leastsq(err,pE[:],args=(xind,yVR))
#    VR = itp.UnivariateSpline(diff_time[good_ind],VnR[good_ind,i],s=0)
#    VnR_sm.append(VR)
##        VnR_sm[:,j,i] = VR
##        yVX = [VnX[good_ind[j],i],VnX[good_ind[j+1],i]]
##        VX,success = opt.leastsq(err,pE[:],args=(xind,yVX))
#    VX = itp.UnivariateSpline(diff_time[good_ind],VnX[good_ind,i],s=0)
#    VnX_sm.append(VX)
##        VnX_sm[:,j,i] = VX
##        yIR = [InR[good_ind[j],i],InR[good_ind[j+1],i]]
##        IR,success = opt.leastsq(err,pE[:],args=(xind,yIR))
#    IR = itp.UnivariateSpline(diff_time[good_ind],InR[good_ind,i],s=0)
#    InR_sm.append(IR)
##        InR_sm[:,j,i] = IR
##        yIX = [InX[good_ind[j],i],InX[good_ind[j+1],i]]
##        IX,success = opt.leastsq(err,pE[:],args=(xind,yIX))
#    IX = itp.UnivariateSpline(diff_time[good_ind],InX[good_ind,i],s=0)
#    InX_sm.append(IX)
##        InX_sm[:,j,i] = IX
##        diff_time_sm[j]=diff_time[good_ind[j]]
##diff_time_sm[-1] = diff_time[good_ind[-1]]

#print shape(InR_sm)
#ant_time = arange(diff_time[0]-1,diff_time[-1]+1,0.01)
#pylab.plot(diff_time[good_ind],Temp_gain[good_ind,1000])
#pylab.plot(ant_time,Temp_gain_sm[1000](ant_time))
#pylab.savefig(data_dir+'cal_temp',dpi=300)
#pylab.clf()


R_ant,X_ant,F_ant = fc.imped_skrf(ant_s11_file,0.0)
R_amp,X_amp,F_amp = fc.imped_skrf(amp_s_file,0.0)
Effic = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
Eff_sm = fc.smooth(Effic,F_ant,0.01)
R_ant_sm = fc.smooth(R_ant,F_ant,0)
X_ant_sm = fc.smooth(X_ant,F_ant,0)
R_amp_sm = fc.smooth(R_amp,F_amp,4e3)
X_amp_sm = fc.smooth(X_amp,F_amp,4e3)
Z_amp = R_amp_sm(new_freq)+1j*X_amp_sm(new_freq)
Z_ant = R_ant_sm(new_freq)+1j*X_ant_sm(new_freq)
Z50 = 50.*ones(len(new_freq))
Z100 = 100.*exp(2*pi*array(new_freq)*400*1e-6*1j)

timescale = 32
binscale = 32
full_rebin_ant = []
full_rebin_time = []
full_rebin_freq = []

data_arrays = os.listdir(data_dir)
#dates = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
dates = sys.argv[1]

selected_ind = []
for f in range(0,len(data_arrays)):
    fname = data_arrays[f]
    date = 'June'+dates
    if fname.split('-')[0]==date:
	if fname.split('_')[-1]=='antenna.txt':
            selected_ind.append(f)
	

# want to test with a single date.
#for date_ind in range(1,2):
for s in range(0,len(selected_ind)):
    fname = data_arrays[selected_ind[s]]
#    fname = 'June'+dates[date_ind]+'_antenna.txt'
    if fname.split('_')[-1]=='antenna.txt':
        print 'Data Collection Date is:', fname.split('_')[0]
	sub_date = fname.split('_')[0]
        rebin_ant = loadtxt(data_dir+fname)
        for fname2 in data_arrays:
            if fname2.split('_')[0]==fname.split('_')[0]:
                if fname2.split('_')[-1]=='mask.txt':
                    ant_mask = loadtxt(data_dir+fname2)
                elif fname2.split('_')[-1]=='time.txt':
                    ant_time = loadtxt(data_dir+fname2)
                elif fname2.split('_')[-1]=='freq.txt':
                    new_freq = loadtxt(data_dir+fname2)
        new_freq = array(new_freq)
        rebin_ant = array(rebin_ant)
        ant_mask= array(ant_mask)
        ant_mask_full = zeros((len(ant_mask),len(ant_mask[0])))
        for j in range(0,len(new_freq)):
            test_mask = fc.timeflag(rebin_ant[:,j],ant_mask[:,j],ant_time,3.,30)
            ant_mask_full[:,j] = test_mask

        data_eff_cal = []
        for k in range(0,len(rebin_ant)):
            lin_data = 10.0**(rebin_ant[k]/10.0)
            corr_data = fc.effcal(Eff_sm,new_freq,lin_data)
            data_eff_cal.append(corr_data)
#        data_eff_cal_db = 10.0*log10(data_eff_cal)
        data_eff_cal_db = rebin_ant

        ant_eff_noise_corr = []
        for i in range(0,len(data_eff_cal_db)):
##            index1 = 0
##            index2 = 0 
##            time_comp = 100.
##            corr_index = 0
##            for j in range(0,len(diff_time)):
##                if diff_time[j]<ant_time[i]:
##                    time_comp = abs(diff_time[j]-ant_time[i])
##                    index2 = index1
##                    index1 = j
##                    corr_index = index1
##                    if abs(diff_time[index2]-ant_time[i])>0.5:
##                        index2 = j
##                elif diff_time[-1]<ant_time[i]:
##                    corr_index = -1
##            if index1==index2:
##                Psky = fc.noise_corr(data_eff_cal_db[i],Vn[index1],In[index1],new_freq,Z_amp,Z_ant,Gain[index1],Temp_gain[index1])
##                Tsky = Temp_gain[index1]*Psky
##                ant_eff_noise_corr.append(Tsky)
##            else:
##                Psky = fc.noise_corr(data_eff_cal_db[i],(Vn[index1]+Vn[index2])/2.,(In[index1]+In[index2])/2.,new_freq,Z_amp,Z_ant,(Gain[index1]+Gain[index2])/2.,(Temp_gain[index1]+Temp_gain[index2])/2.)
##                Tsky = (Temp_gain[index1]+Temp_gain[index2])*Psky/2.
##                ant_eff_noise_corr.append(Tsky)
##            if index2>index1:
##                corr_index = index1
##            elif index1>index2:
##                corr_index = index2
##            elif index1==index2:
##                corr_index = index1
            Tsky = zeros(len(data_eff_cal_db[0]))
            for j in range(0,len(data_eff_cal_db[0])):
#               TGain_s = Temp_gain_sm[j](ant_time[i])
#               Vn_s = VnR_sm[j](ant_time[i])+1j*VnX_sm[j](ant_time[i])
#               In_s = InR_sm[j](ant_time[i])+1j*InX_sm[j](ant_time[i])
#               Gain_s = GainR_sm[j](ant_time[i])+1j*GainX_sm[j](ant_time[i])
##               TGain_s= Temp_gain_sm[0,corr_index,j]*ant_time[i]+Temp_gain_sm[1,corr_index,j]
##               Vn_s = VnR_sm[0,corr_index,j]*ant_time[i]+VnR_sm[1,corr_index,j]+VnX_sm[0,corr_index,j]*ant_time[i]*1j+VnX_sm[1,corr_index,j]*1j
##               In_s = InR_sm[0,corr_index,j]*ant_time[i]+InR_sm[1,corr_index,j]+InX_sm[0,corr_index,j]*ant_time[i]*1j+InX_sm[1,corr_index,j]*1j
##               Gain_s = GainR_sm[0,corr_index,j]*ant_time[i]+GainR_sm[1,corr_index,j]+GainX_sm[0,corr_index,j]*ant_time[i]*1j+GainX_sm[1,corr_index,j]*1j
               TGain_s = Temp_gain[j]
               Vn_s = VnR[j]+1j*VnX[j]
               In_s = InR[j]+1j*InX[j]
               Gain_s = GainR[j]+1j*GainX[j]
               Tsky[j] = fc.noise_corr(data_eff_cal_db[i,j],Vn_s,In_s,new_freq[j],Z_amp[j],
                                    Z_ant[j],Gain_s,TGain_s)
#               Tsky[j] = TGain_s*Psky
            ant_eff_noise_corr.append(Tsky)
                

        ant_eff_noise_corr = array(ant_eff_noise_corr)
#        median_time_eff = ma.median(ant_eff_noise_corr,axis=1)
        single_time_eff = ant_eff_noise_corr[:,4000]
	full_single_median = ma.median(single_time_eff)

        single_time_lim = []
        for i in range(0,len(single_time_eff)):
            if single_time_eff[i]>full_single_median/5.:
		if single_time_eff[i]<5.*full_single_median:
                    single_time_lim.append(single_time_eff[i])            
        median_single_time = ma.median(single_time_lim)
        std_single_time = ma.std(single_time_lim)
#        total_median = ma.median(median_time_eff)
	print 'Plotting Time Variation (median is %0.1f)...' %median_single_time 
	
#	fig1 = plt.figure()
        pylab.scatter(ant_time,single_time_eff,c='b',edgecolor='b',s=3)
        pylab.scatter(ant_time,ones(len(ant_time))*(median_single_time),c='g',edgecolor='g',s=3)
        pylab.scatter(ant_time,ones(len(ant_time))*(median_single_time+3*std_single_time),c='r',edgecolor='r',s=3)
        pylab.scatter(ant_time,ones(len(ant_time))*(median_single_time-3*std_single_time),c='c',edgecolor='c',s=3)
	pylab.legend(('Data','Median','+3 Sigma','-3 Sigma'))
	pylab.grid()
#	pylab.ylim(-2e4,2e4)
	ymin,ymax = pylab.ylim()
	if ymin<0:
	    pylab.ylim(ymin=0)
	if ymax>5e3:
	    pylab.ylim(ymax=5e3)
        pylab.xlim(ant_time[0],ant_time[-1]+0.25)
	pylab.xlabel('Time (Hours since Midnight June 01 GMT)')
	pylab.ylabel('Temperature (Kelvin)')
	pylab.title('Time Variability of Data at %0.1f MHz' %new_freq[4000])
	pylab.savefig(out_dir+sub_date+'_single_freq_median_avgcal',dpi=300)
#        pylab.savefig(data_dir+fname.split('_')[0]+'single_freq_median_avgcal',dpi=300)
        pylab.clf()
        
#        print total_median

        new_mask_full = zeros((len(ant_mask_full),len(ant_mask_full[0])))
        for i in range(0,len(ant_mask_full)):
            nandata = where(isnan(ant_eff_noise_corr[i]))[0]
            nandata = array(nandata)
            for k in range(0,len(nandata)):
                index = nandata[k]
                new_mask_full[i,index]=1.0
                ant_eff_noise_corr[i,index] = 0.0
            new_mask_single = fc.flagging(ant_eff_noise_corr[i],new_freq,3.,30)
            new_mask_full[i] = new_mask_single
#            spike_mask_single = fc.spike_flag(ant_eff_noise_corr[i],100)
            for j in range(0,len(ant_mask_full[0])):
                if ant_mask_full[i,j]==1.0:
                    new_mask_full[i,j] = 1.0
                if single_time_eff[i]>(median_single_time+3*std_single_time):
                    new_mask_full[i,j] = 1.0
                elif single_time_eff[i]<(median_single_time-3*std_single_time):
                    new_mask_full[i,j] = 1.0
                elif single_time_eff[i]<0:
                    new_mask_full[i,j] = 1.0
#                if single_time_eff[i]<=1000.:
#                    new_mask_full[i,j] = 1.0
#                if spike_mask_single[j]==1.0:
#                    new_mask_full[i,j] = 1.0
#	    nandata = where(isnan(ant_eff_noise_corr[i]))
#	    nandata = array(nandata)
#	    for k in range(0,len(nandata[0])):
#                index = nandata[0,i]
#		new_mask_full[i,index]=1.0
                
	
	print 'Plotting Mask...'
        pylab.imshow(new_mask_full,vmin=-0.1,vmax=1.1,aspect=100./len(new_mask_full),extent=(40,140,len(new_mask_full),0.))
        pylab.ylabel('Time (Sample)')
        pylab.xlabel('Freq (MHz)')
        pylab.colorbar()
        pylab.title('Flagged Frequencies and Times for Data')
	pylab.savefig(out_dir+sub_date+'_data_mask',dpi=300)
#        pylab.savefig(data_dir+fname.split('_')[0]+'_data_mask',dpi=300)
        pylab.clf()
                    
	print 'Plotting Single Time Masking...'
        single_time = int(len(new_mask_full)/2)
	time_used = single_time
	okay = 0
#	for i in range(0,single_time-3):
#	    if okay==0:
#		if ma.mean(new_mask_full[single_time+i])<0.5:
#		    okay=1
#		    time_used = single_time+i	
	for j in range(0,single_time):
	    if okay==0:
		print 'Bad File is:',time_used
		meanmask = ma.median(new_mask_full[time_used])
		time_used = time_used-j
		if meanmask<0.9:
		    okay=1
#		    time_used = single_time-j
	pylab.scatter(new_freq,ant_eff_noise_corr[time_used],c='b',edgecolor='b',s=3)
	masked_single_time = ma.array(ant_eff_noise_corr[time_used],mask=new_mask_full[time_used])
	compressed_single_time = ma.compressed(masked_single_time)
	masked_single_freq = ma.array(new_freq,mask=new_mask_full[time_used])
	compressed_single_freq = ma.compressed(masked_single_freq)
	print 'Compressed Shape is:',len(compressed_single_time)
	pylab.scatter(compressed_single_freq,compressed_single_time,c='g',edgecolor='g',s=3)
	pylab.grid()
	pylab.xlabel('Frequency (MHz)')
	pylab.xlim(40,140)
	pylab.ylabel('Temperature (Kelvin)')
	pylab.legend(('Unmasked Sample','Masked Sample'))
	pylab.title('Example Masking Dataset')
#	pylab.ylim(-5e4,5e4)
	pylab.ylim(0,2e4)
	pylab.savefig(out_dir+sub_date+'_single_time_masking',dpi=300)
	pylab.clf()


        freqbin_full = []
        freqbin_mask_full = []
        for i in range(0,len(ant_eff_noise_corr)):
            freqbin_data,freqbin_mask,freqbin_freq = fc.rebin(ant_eff_noise_corr[i],
                                                              new_mask_full[i],
                                                              new_freq,binscale)
            freqbin_full.append(freqbin_data)
            freqbin_mask_full.append(freqbin_mask)
        minimum_time = ant_time[0]
        time_ind = 0
	min_time_ind = 0
        timebin_time = []
	print 'Shape of Time Data is:',shape(ant_time)
        for i in range(0,len(ant_eff_noise_corr)):
	    if time_ind<timescale:
		if time_ind<1:
		    minimum_time = ant_time[i]
		    time_ind = time_ind+1
                elif absolute(ant_time[i]-minimum_time)<1./12:
		    minimum_time = ant_time[i]
                    time_ind = time_ind + 1
		else:
                    print 'Number of Time Samples to Bin:',time_ind
                    timebin_data,timebin_mask = fc.timerebin(freqbin_full[min_time_ind:i],
                                                             freqbin_mask_full[min_time_ind:i])
                    full_rebin_time.append(ma.mean(ant_time[min_time_ind:i]))
                    full_rebin_ant.append(timebin_data)
                    min_time_ind = i+1
                    minimum_time = ant_time[i]
                    time_ind = 0
	    else:
		print 'Number of Time Samples to Bin:',time_ind
                timebin_data,timebin_mask = fc.timerebin(freqbin_full[min_time_ind:i],
                                                         freqbin_mask_full[min_time_ind:i])
                full_rebin_time.append(ma.mean(ant_time[min_time_ind:i]))
                full_rebin_ant.append(timebin_data)
                min_time_ind = i+1
                minimum_time = ant_time[i]
                time_ind = 0
        if time_ind>=5:
#        if (len(ant_eff_noise_corr)-time_ind)>5:
            timebin_data,timebin_mask = fc.timerebin(freqbin_full[min_time_ind:-1],
                                                     freqbin_mask_full[min_time_ind:-1])
            full_rebin_time.append(ma.mean(ant_time[min_time_ind:-1]))
            full_rebin_ant.append(timebin_data)
        print 'Number of Excess Files:',time_ind
        full_rebin_freq = freqbin_freq
        print 'New Data Array Shape:',shape(full_rebin_ant)

full_rebin_ant = array(full_rebin_ant)
for i in range(0,len(full_rebin_ant)):
    for j in range(0,len(full_rebin_ant[0])):
     	if isnan(full_rebin_ant[i,j])==True:
            full_rebin_ant[i,j]=0.

full_rebin_time = array(full_rebin_time)
full_rebin_freq = array(full_rebin_freq)
savetxt(out_dir+date+'_processed_data_avgcal.txt',full_rebin_ant,delimiter = ' ')
savetxt(out_dir+date+'_processed_time_avgcal.txt',full_rebin_time,delimiter = ' ')
savetxt(out_dir+date+'_processed_freq_avgcal.txt',full_rebin_freq,delimiter = ' ')
full_rebin_ant = []
full_rebin_time = []
full_rebin_freq = []
                
        

        

        
