from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import data_analysis_funcs as fc

#directory = 'converted_VNA_data/'
#ant_data = loadtxt(directory+'log_mag_good_cable.txt')
#ant_data = loadtxt(directory+'long_ant.txt')
#ant_data = loadtxt(directory+'short_ant.txt')
#amp_data = loadtxt(directory+'m_amp1_12v.txt')
#amp_data = loadtxt(directory+'amp_alg.txt')
#R_ant, X_ant, F_ant = fc.imped(ant_data,0.0)
#R_amp, X_amp, F_amp = fc.imped(amp_data,0.0)
#E = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
#E_sm = fc.smooth(E,F_ant,0.3)


#pylab.plot(F_ant,E)
#pylab.savefig('efficiency_test_phase_5_2.jpg',dpi=200)
#pylab.clf()

directory2 = 'Isla_Guadalupe_data/long_data_village/'
#directory2 = 'Algonquin_data/raw_data/log_period_amp_at_zenith/'
#directory2 = 'Algonquin_data/raw_data/long_test/'
dirList = os.listdir(directory2)
load_data = []
load_time = []
short_data = []
short_time = []
term_data = []
term_time = []
noise_data = []
noise_time = []
load_volt = []
short_volt = []
term_volt = []
noise_volt = []

#i = 0
for fname in dirList:
    test_file = directory2 + fname
#    if i>10:
    if fname!='plots':
        if fname.split('_')[3]=='50ohm.dat':
            ltime,form,ldata,lmask,lfreq = fc.loadsingle(test_file)
            load_time.append(ltime)
            load_data.append(ldata)
            f = open(test_file)
            test = f.readlines(1)
            for line in range(0,len(test)):
                if test[line].split('\t')[0]=='#  VIN:':
                    lim = test[line].split('\t')[1]
                    load_volt.append(float(lim.split('\n')[0]))
        elif fname.split('_')[3]=='short.dat':
            ltime,form,ldata,lmask,lfreq = fc.loadsingle(test_file)
            short_time.append(ltime)
            short_data.append(ldata)
            f = open(test_file)
            test = f.readlines(1)
            for line in range(0,len(test)):
                if test[line].split('\t')[0]=='#  VIN:':
                    lim = test[line].split('\t')[1]
                    short_volt.append(float(lim.split('\n')[0]))
        elif fname.split('_')[3]=='noise.dat':
            ltime,form,ldata,lmask,lfreq = fc.loadsingle(test_file)
            noise_data.append(ldata)
            noise_time.append(ltime)
            f = open(test_file)
            test = f.readlines(1)
            for line in range(0,len(test)):
                if test[line].split('\t')[0]=='#  VIN:':
                    lim = test[line].split('\t')[1]
                    noise_volt.append(float(lim.split('\n')[0]))
        elif fname.split('_')[3]=='open.dat':
            ltime,form,ldata,lmask,lfreq = fc.loadsingle(test_file)
            term_data.append(ldata)
            term_time.append(ltime)
            f = open(test_file)
            test = f.readlines(1)
            for line in range(0,len(test)):
                if test[line].split('\t')[0]=='#  VIN:':
                    lim = test[line].split('\t')[1]
                    term_volt.append(float(lim.split('\n')[0]))
        
        
#    i+=1

width = 250.0/len(ldata)
freq_data = arange(0,250,width)
load_data = array(load_data)
short_data = array(short_data)
noise_data = array(noise_data)
term_data = array(term_data)
term_volt = array(term_volt)
noise_volt = array(noise_volt)
short_volt = array(short_volt)
load_volt = array(load_volt)
l_sort = argsort(load_volt)
s_sort = argsort(short_volt)
n_sort = argsort(noise_volt)
t_sort = argsort(term_volt)
ld_sort = zeros((len(load_data),len(load_data[0])))
sd_sort = zeros((len(short_data),len(short_data[0])))
nd_sort = zeros((len(noise_data),len(noise_data[0])))
td_sort = zeros((len(term_data),len(term_data[0])))
for i in range(0,len(l_sort)):
    ld_sort[i,:] = load_data[l_sort[i],:]
for i in range(0,len(s_sort)):
    sd_sort[i,:] = short_data[s_sort[i],:]
for i in range(0,len(n_sort)):
    nd_sort[i,:] = noise_data[n_sort[i],:]
for i in range(0,len(t_sort)):
    td_sort[i,:] = term_data[t_sort[i],:]
lv_sort = sort(load_volt)
tv_sort = sort(term_volt)
nv_sort = sort(noise_volt)
sv_sort = sort(short_volt)

fitfunc = lambda p,x: p[0]+p[1]*x+p[2]*x*x
errfunc = lambda p,x,y: fitfunc(p,x)-y
pl = zeros((3,len(load_data[0])))
pn = zeros((3,len(noise_data[0])))
ps = zeros((3,len(short_data[0])))
pt = zeros((3,len(term_data[0])))
p0 = [-145.,11.1,-0.5]


#E_fifty = fc.effic([50.0,],[0.0,],freq_data,R_amp,X_amp,F_amp)
#E_fsm = fc.smooth(E_fifty,F_amp,0.475)

#Gain_avgc, Gain_c = fc.calcgain(load_data,short_data,E_fsm(freq_data))
#G_avgc_sm = fc.smooth(Gain_avgc,freq_data,3e-20)
#short_time = array(short_time)
#short_data = array(short_data)

#pylab.plot(freq_data,Gain_avgc)
#pylab.savefig('guad_prelim_gain.jpg',dpi=200)
#pylab.clf()

#pylab.plot(freq_data,E_sm(freq_data))
#pylab.ylim(-1e-13,1e-13)
#pylab.savefig('gain_smoothed.jpg',dpi=200)
#pylab.clf()

fnum = 0
full_data = []
mask_full = []
time_full = []
uneff_data = []
lim_data = []
lim_mask = []
lim_time = []
freqscale = 60
timescale = 30


#for fname in dirList:
#    if fname.split('_')[3]=='antenna.dat':
#        test_file = directory2+fname
#        time, form, sub_data,mask_data,freq_data = fc.loadsingle(test_file)
#        sub_data = 10.0**(sub_data/10.0)
#        width = 250.0/len(sub_data)
#        freq_data = arange(0,250,width)
#        mask_data = fc.flagging(sub_data,freq_data,3.0)
#        diff = 24.0
#        for j in range(0,len(short_time)):
#            if abs(short_time[j]-time)<diff:
#                diff= abs(short_time[j]-time)
#                index = j
#        sh_used = short_data[index,:]
#        sh_used = 10.0**(sh_used/10.0)
#        corr_data = fc.gaincal(G_avgc_sm,freq_data,sub_data,sh_used)
#        uneff_data.append(corr_data)
#        corr_data = fc.effcal(E_sm,freq_data,corr_data)
#        new_data, new_mask, new_freq = fc.rebin(sub_data,mask_data,freq_data,freqscale)
#        full_data.append(sub_data)
#        time_full.append(time)
#        mask_full.append(mask_data)
#        minfreq = where(new_freq<=30.0)[0][-1]
#        maxfreq = where(new_freq>=249.0)[0][0]
#        lim_data.append(new_data[minfreq:maxfreq])
#        lim_mask.append(new_mask[minfreq:maxfreq])
#        lim_time.append(time)
#        print fnum
#        fnum+=1
#        if fnum==timescale:
#            lim_data = array(lim_data)
#            lim_mask = array(lim_mask)
#            lim_time = array(lim_time)
#            nl_data, nl_mask, nl_time = fc.timerebin(lim_data,lim_mask,lim_time)
#            lim_data = []
#            lim_mask = []
#            lim_time = []
#            full_data.append(nl_data)
#            mask_full.append(nl_mask)
#            time_full.append(nl_time)
#            print nl_time, nl_data[len(new_freq[minfreq:maxfreq])/2]
#            fnum = 0
#        print minfreq,maxfreq
#        p_fit = fc.modelfit(corr_data[minfreq:maxfreq],
#                            freq_data[minfreq:maxfreq],'lin_cos',
#                            mask_data[minfreq:maxfreq])
#        pylab.scatter(freq_data,corr_data,s=5,edgecolor='b')
#        fitfunc = lambda p,x: p[0]*cos(2*pi*p[1]*x+p[2]%(2*pi))+ p[3] #cosine + offset        
#        fitfunc = lambda p,x: 10**(p[0]*log10(x)+p[1])*(cos(2*pi*p[2]*x+p[3]%(2*pi))+p[4])#linear + sinusoid
#        p_0 = [-1.25,-5.7,0.051,1.9,2.4]
#        fit_data = fitfunc(p_fit,freq_data)
#        print fit_data[5000],corr_data[5000]
#        pylab.plot(freq_data,fitfunc(p_0,freq_data),'r-o',ms=3,markeredgecolor='r')
#        pylab.plot(freq_data,fit_data,'g-o',ms=3,markeredgecolor='g')
#        pylab.xlim(30.0,100.0)
#        pylab.xlabel('Frequency (MHz)')
#        pylab.ylim(0.0,1e-7)
#        pylab.ylabel('Temperature (K)')
#        pylab.savefig('fit_log_period_cal_'+str(fnum)+'.jpg',dpi=200)
#        pylab.clf()
    
#full_data = array(full_data)
#print shape(full_data)
#full_data = 10*log10(full_data)
#mask_full = array(mask_full)
#freq_full = new_freq[minfreq:maxfreq]
#freq_full = array(freq_full)
#time_full = array(time_full)
#extra_flag = zeros((len(mask_full),len(mask_full[0])))

#for f in range(0,len(freq_full)):
#    extra_flag[:,f] = fc.timeflag(full_data[:,f],mask_full[:,f],time_full,
#                                  freq_full[f],3.0)
#    ratio = sum(extra_flag[:,f])/len(extra_flag)
#    print freq_full[f], ratio

#comp_data = ma.array(full_data,mask=extra_flag)

#fc.waterfallplot(full_data,[-90,-50],time_full,freq_full,
#                 'Uncalibrated Data','guad_test_long.jpg')
#fc.waterfallplot(10**(full_data/10),[0,1e-6],time_full,freq_full,
#                 'Uncalibrated Data','guad_test_long_lin.jpg')

#mask_mean = sum(mask_full,axis=0)/len(mask_full)
#mean_data = ma.mean(full_data,axis=0)
#minfreq = where(freq_data<=30.0)[0][-1]
#maxfreq = where(freq_data>=100.0)[0][0]
#p_fit = fc.modelfit(mean_data[minfreq:maxfreq],
#                    freq_data[minfreq:maxfreq],'lin_cos',
#                    mask_full[0,minfreq:maxfreq])
#fitfunc = lambda p,x: 10**(p[0]*log10(x)+p[1])*(cos(2*pi*p[2]*x+p[3]%(2*pi))+p[4])#linear + sinusoid
#p_0 = [-1.25,-5.7,0.051,1.9,2.4]
#fit_data = fitfunc(p_fit,freq_data)
#all_mask = where(mask_mean>0.0)[0]
#limmask = where(abs(mean_data-fitfunc(p_0,freq_data))>2e-8)[0]
#newmask = zeros(len(mean_data))
#newmask[all_mask] = 1.0
#newmask[limmask] = 1.0
#newmask[0:minfreq] = 1.0
#newmask[maxfreq:-1] = 1.0
#new_data = ma.array(mean_data,mask=newmask)
#new_data = new_data.compressed()
#new_freq = ma.array(freq_data,mask=newmask)
#new_freq = new_freq.compressed()
#p_new = fc.modelfit(new_data,new_freq,'lin_cos',zeros(len(new_data)))
#pylab.scatter(freq_data,mean_data,s=5,edgecolor='b')
#pylab.plot(freq_data,fitfunc(p_0,freq_data),'r-o',ms=1,markeredgecolor='r')
#pylab.plot(freq_data,fit_data,'g-o',ms=1,markeredgecolor='g')
#pylab.plot(new_freq,fitfunc(p_new,new_freq),'c-o',ms=1,markeredgecolor='c')
#pylab.xlim(30.0,100.0)
#pylab.xlabel('Frequency (MHz)')
#pylab.ylim(0.0,1e-7)
#pylab.savefig('fit_log_period_cal_mean.jpg',dpi=200)
#pylab.clf()

#pylab.plot(freq_data, mask_mean)
#pylab.plot(freq_data,newmask,c='g')
#pylab.ylim(-0.1,1.1)
#pylab.savefig('mask_mean.jpg',dpi=200)
#pylab.clf()

#pylab.plot(freq_data,mean_data-fitfunc(p_fit,freq_data),'g')
#pylab.plot(freq_data,mean_data-fitfunc(p_0,freq_data),'r')
#pylab.plot(freq_data,mean_data-fitfunc(p_new,freq_data),'c')
#pylab.xlim(30.0,100.0)
#pylab.xlabel('Frequency (MHz)')
#pylab.ylim(-2e-8,2e-8)
#pylab.savefig('mean_fit_residuals.jpg',dpi=200)
#pylab.clf()
           
#amp_header = []
#Rs = []
#Xs = []
#Fs = []
#Es = []
#ant_data = loadtxt(directory+'short_ant.txt')
#ant_data = loadtxt(directory+'med_ant.txt')
#ant_data = loadtxt(directory+'long_ant.txt')
#R_ant,X_ant,F_ant = fc.imped(ant_data,0.0)
#pylab.plot(F_ant, R_ant)
#pylab.xlim(0.0,250.0)
#pylab.savefig('short_ant_real_impedence.jpg',dpi=200)
#pylab.clf()
#pylab.plot(F_ant,X_ant)
#pylab.xlim(0.0,250.0)
#pylab.savefig('short_ant_imag_impedence.jpg',dpi=200)
#pylab.clf()

#for fname in dirList:
#    test_file = directory+fname
#    time, form, sub_data,mask_data,freq_data = fc.loadsingle(test_file)

#    width = 125.0/len(sub_data)
#    freq_data = arange(0,125,width)

#    sub_data = 10.0**(sub_data/10.0)

#    mask_data = fc.flagging(sub_data,freq_data,3.0)

#    print ma.median(sub_data), sum(mask_data)
#    if fname.split('_')[0]=='amp':
#        amp_name = fname.split('.')[0]
#        amp_header.append(amp_name)
#        amp_data=loadtxt(directory+fname)
#        R_amp, X_amp, F_amp = fc.imped(amp_data,0.0)
#        E = fc.effic(R_ant,X_ant,F_ant,R_amp,X_amp,F_amp)
#        Rs.append(R_amp)
#        Xs.append(X_amp)
#        Fs.append(F_amp)
#        Es.append(E)

#for r in range(0,7):
#    pylab.plot(Fs[r],Rs[r],label=str(amp_header[r]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(0,500)
#pylab.savefig('amplifiers_real_impedence.jpg',dpi=200)
#pylab.clf()

#for r in range(7,14):
#    pylab.plot(Fs[r],Rs[r],label=str(amp_header[r]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(0,500)
#pylab.savefig('amplifiers_real_impedence_2.jpg',dpi=200)
#pylab.clf()

#for r in range(14,len(Rs)):
#    pylab.plot(Fs[r],Rs[r],label=str(amp_header[r]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(0,500)
#pylab.savefig('amplifiers_real_impedence_3.jpg',dpi=200)
#pylab.clf()

#for x in range(0,7):
#    pylab.plot(Fs[x],Xs[x],label=str(amp_header[x]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(-500,500)
#pylab.savefig('amplifiers_imag_impedence.jpg',dpi=200)
#pylab.clf()
        
#for x in range(7,14):
#    pylab.plot(Fs[x],Xs[x],label=str(amp_header[x]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(-500,500)
#pylab.savefig('amplifiers_imag_impedence2.jpg',dpi=200)
#pylab.clf()
        
#for x in range(14,len(Xs)):
#    pylab.plot(Fs[x],Xs[x],label=str(amp_header[x]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(-500,500)
#pylab.savefig('amplifiers_imag_impedence3.jpg',dpi=200)
#pylab.clf()

#for e in range(0,7):
#    pylab.plot(Fs[e],Es[e],label=str(amp_header[e]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(0,1.0)
#pylab.savefig('amplifiers_eff_short.jpg',dpi=200)
#pylab.clf()
        
#for e in range(7,14):
#    pylab.plot(Fs[e],Es[e],label=str(amp_header[e]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(0,1.0)
#pylab.savefig('amplifiers_eff_short2.jpg',dpi=200)
#pylab.clf()
        
#for e in range(14,len(Es)):
#    pylab.plot(Fs[e],Es[e],label=str(amp_header[e]),linestyle='dashed',linewidth=8)
#pylab.legend()
#pylab.xlim(0.0,350.0)
#pylab.ylim(0,1.0)
#pylab.savefig('amplifiers_eff_short3.jpg',dpi=200)
#pylab.clf()
        
