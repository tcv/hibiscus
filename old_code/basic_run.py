from numpy import *
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
from scipy import optimize
import os
import data_analysis_funcs as fc


directory2 = 'Isla_Guadalupe_data/long_data_village/'
#directory2 = 'Algonquin_data/raw_data/log_period_amp_at_zenith/'
#directory2 = 'Algonquin_data/raw_data/long_test/'
dirList = os.listdir(directory2)
volt_data = loadtxt('voltage_values_for_long_data.txt',delimiter = ' ')
volt_data = array(volt_data)
fnum = 0
full_data = []
mask_full = []
time_full = []
volt_full = []
uneff_data = []
lim_data = []
lim_mask = []
lim_time = []
freqscale = 50
timescale = 1
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

for fname in dirList:
    test_file = directory2 + fname
    if fname!='plots':
        if fname.split('_')[3]=='50ohm.dat':
            ltime,form,ldata,lmask,lfreq,lvolt = fc.loadsingle(test_file)
            load_time.append(ltime)
            load_data.append(ldata)
            load_volt.append(lvolt)
        elif fname.split('_')[3]=='short.dat':
            ltime,form,ldata,lmask,lfreq,lvolt = fc.loadsingle(test_file)
            short_time.append(ltime)
            short_data.append(ldata)
            short_volt.append(lvolt)
        elif fname.split('_')[3]=='noise.dat':
            ltime,form,ldata,lmask,lfreq,lvolt = fc.loadsingle(test_file)
            noise_data.append(ldata)
            noise_time.append(ltime)
            noise_volt.append(lvolt)
        elif fname.split('_')[3]=='open.dat':
            ltime,form,ldata,lmask,lfreq,lvolt = fc.loadsingle(test_file)
            term_data.append(ldata)
            term_time.append(ltime)
            term_volt.append(lvolt)

width = 250.0/len(ldata)
freq_data = arange(0,250,width)

load_data = array(load_data)
short_data = array(short_data)
#noise_data = array(noise_data)
#term_data = array(term_data)
load_time = array(load_time)
short_time = array(short_time)

#term_volt = array(term_volt)
noise_volt = array(noise_volt)
short_volt = array(short_volt)
#load_volt = array(load_volt)

eff = ones(len(load_data[0]))
load_mod = load_data
short_mod = short_data
volt_l =zeros(len(load_data))
volt_s = zeros(len(short_data))
for j in range(0,len(load_data)):
    volt_l[j] = volt_data[1,0]
    volt_s[j] = volt_data[1,0]
    for i in range(0,len(volt_data[0])):
        if volt_data[0,i]<=load_time[j]:
            volt_l[j] = volt_data[1,i]
            if volt_l[j]<0.0:
                volt_l[j] = 0.0
        if volt_data[0,i]<=short_time[j]:
            volt_s[j] = volt_data[1,i]
            if volt_s[j]<0.0:
                volt_s[j] =0.0
    print load_time[j], volt_l[j]
                
#l_sort = argsort(load_volt)
#s_sort = argsort(short_volt)
#n_sort = argsort(noise_volt)
#t_sort = argsort(term_volt)

#ld_sort = zeros((len(load_data),len(load_data[0])))
#sd_sort = zeros((len(short_data),len(short_data[0])))
#nd_sort = zeros((len(noise_data),len(noise_data[0])))
#td_sort = zeros((len(term_data),len(term_data[0])))

#for i in range(0,len(l_sort)):
#    ld_sort[i,:] = load_data[l_sort[i],:]
#for i in range(0,len(s_sort)):
#    sd_sort[i,:] = short_data[s_sort[i],:]
#for i in range(0,len(n_sort)):
#    nd_sort[i,:] = noise_data[n_sort[i],:]
#for i in range(0,len(t_sort)):
#    td_sort[i,:] = term_data[t_sort[i],:]

#lv_sort = sort(load_volt)
#tv_sort = sort(term_volt)
#nv_sort = sort(noise_volt)
#sv_sort = sort(short_volt)

fitfunc = lambda p,x: p[0]+p[1]*x+p[2]*x*x
errfunc = lambda p,x,y: fitfunc(p,x)-y
#pl = zeros((3,len(load_data[0])))
#pn = zeros((3,len(noise_data[0])))
#ps = zeros((3,len(short_data[0])))
#pt = zeros((3,len(term_data[0])))
#p0 = [-145.,11.1,-0.5]

#for i in range(0,len(load_data[0])):
#    pl[:,i], success = optimize.leastsq(errfunc,p0[:],args=(lv_sort,ld_sort[:,i]))
#    pt[:,i], success = optimize.leastsq(errfunc,p0[:],args=(tv_sort,td_sort[:,i]))
#    pn[:,i], success = optimize.leastsq(errfunc,p0[:],args=(nv_sort,nd_sort[:,i]))
#    ps[:,i], success = optimize.leastsq(errfunc,p0[:],args=(sv_sort,sd_sort[:,i]))
#    print i, pl[1,i],pt[1,i],pn[1,i],ps[1,i]

p0 = [-145.0,12.0,-0.45]
for i in range(0,len(load_mod)):
    load_mod[i,:] = load_mod[i,:]*fitfunc(p0,11.25)/fitfunc(p0,volt_l[i])
    short_mod[i,:] = short_mod[i,:]*fitfunc(p0,11.25)/fitfunc(p0,volt_s[i])

Gavg_con, G_con = fc.calcgain(load_mod,short_mod,eff)
print shape(G_con)
G_con = array(G_con)
G_con_sm = []
for i in range(0,len(G_con)):
    G_con_sm.append(fc.smooth(G_con[i,:],freq_data,1e-18))
G_con_sm = array(G_con_sm)

print shape(G_con_sm)
#for i in range(0,len(G_con_sm)):
#    pylab.plot(freq_data,G_con_sm[i](freq_data))
#pylab.ylim(0,1e-10)
#pylab.savefig('gain_smoothed_check.jpg',dpi=200)
#pylab.clf()

for fname in dirList:
    if fname!='plots':
        if fname.split('_')[3]=='antenna.dat':
            test_file = directory2+fname
            time, form, sub_data,mask_data,freq_data,volt = fc.loadsingle(test_file)
            for i in range(0,len(volt_data[0])):
                if volt_data[0,i]<=time:
                    volt = volt_data[1,i]
                    if volt<0.0:
                        volt = 0.0
#            print time,volt
            width = 250.0/len(sub_data)
            freq_data = arange(0,250,width)
#            print mean(sub_data)
            p0 = [-145.0,12.0,-0.45]
            for i in range(0,len(freq_data)):
                sub_data[i] = sub_data[i]*fitfunc(p0,11.25)/fitfunc(p0,volt)
#            print mean(sub_data)
            mask_data = zeros(len(sub_data))
#            mask_data = fc.flagging(sub_data,freq_data,3.0)
            sub_data = 10.0**(sub_data/10.0)
            diff = 500.0
            for j in range(0,len(short_time)):
                if abs(short_time[j]-time)<diff:
                    diff = abs(short_time[j]-time)
                    index = j
            print diff
            sh_used = short_mod[index,:]
            sh_used = 10.0**(sh_used/10.0)
            corr_data = fc.gaincal(G_con_sm[index],freq_data,sub_data,sh_used)
            print mean(corr_data)
            new_data, new_mask, new_freq = fc.rebin(corr_data,mask_data,freq_data,freqscale)
#            for i in range(0,len(new_freq)-1):
#                ps_scale = ma.mean(ps[:,i*freqscale:(i+1)*freqscale],axis=1)
#                new_data[i] = new_data[i]*fitfunc(ps_scale,11.25)/fitfunc(ps_scale,volt)
#            ps_scale = ma.mean(ps[:,i*freqscale:-1],axis=1)
#            new_data[-1] = new_data[-1]*fitfunc(ps_scale,11.25)/fitfunc(ps_scale,volt)
#            volt_full.append(volt)
            time_full.append(time)
#        mask_full.append(mask_data)
            minfreq = where(new_freq<=30.0)[0][-1]
            maxfreq = where(new_freq>=249.0)[0][0]
            full_data.append(new_data[minfreq:maxfreq])
#            lim_data.append(new_data[minfreq:maxfreq])
#            lim_mask.append(new_mask[minfreq:maxfreq])
#            lim_time.append(time)
            print fnum
            fnum+=1
#            if fnum==timescale:
#                lim_data = array(lim_data)
#                lim_mask = array(lim_mask)
#                lim_time = array(lim_time)
#                nl_data, nl_mask, nl_time = fc.timerebin(lim_data,lim_mask,lim_time)
#                lim_data = []
#                lim_mask = []
#                lim_time = []
#                full_data.append(nl_data)
#            mask_full.append(nl_mask)
#                time_full.append(nl_time)
#                print (nl_time%24.0), nl_data[len(new_freq[minfreq:maxfreq])/2]
#                fnum = 0
    
full_data = array(full_data)
print shape(full_data)
#full_data = 10*log10(full_data)
#mask_full = array(mask_full)
freq_full = new_freq[minfreq:maxfreq]
freq_full = array(freq_full)
time_full = array(time_full)
#extra_flag = zeros((len(mask_full),len(mask_full[0])))

#for f in range(0,len(freq_full)):
#    extra_flag[:,f] = fc.timeflag(full_data[:,f],mask_full[:,f],time_full,
#                                  freq_full[f],3.0)
#    ratio = sum(extra_flag[:,f])/len(extra_flag)
#    print freq_full[f], ratio

#comp_data = ma.array(full_data,mask=extra_flag)

fc.waterfallplot(full_data,[0,10000],time_full,freq_full,
                 'Calibrated Data','guad_long_flat_volt_adj_gcal.jpg')
#fc.waterfallplot(10**(full_data/10),[0,1e-6],time_full,freq_full,
#                 'Calibrated Data','guad_long_flat_volt_adj_gcal_lin.jpg')

