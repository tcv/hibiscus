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

maindir = '/lustre/tcv/mean_cal_data_old/'
anat_resid = loadtxt('/home/tcv/guad_extras/combined_residuals.dat')
th1 = loadtxt('/home/tcv/guad_extras/theory_one.dat')
th2 = loadtxt('/home/tcv/guad_extras/theory_two.dat')
th3 = loadtxt('/home/tcv/guad_extras/theory_three.dat')
full = False
base_lim = False
gm_lim = False
gu_lim = False
all_type = True
#full dataset:
#if full:
days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14']
two_art = [False,False,True,True,False,False,True,False,False,False,True,True,True,True]
count_a = [4196,90,3340,2506,1888,5444,341,1738,3105,1340,3182,2982,351,2684]
count_b = [0,0,3172,2325,0,0,133,0,0,0,1760,1504,233,1266]
minfgu = [60.,70.,60.,60.,60.,60.,62.,60.,70.,62.,60.,60.,60.,60.]
maxfgu = [90.,90.,90.,90.,90.,90.,77.,90.,75.,80.,80.,90.,80.,90.]
minfgm = [60.,70.,60.,60.,65.,62.,62.,62.,70.,70.,70.,62.,62.,60.]
maxfgm = [90.,75.,75.,75.,70.,77.,77.,77.,75.,70.,75.,85.,77.,75.]
minft = [60.,70.,60.,60.,65.,62.,62.,62.,70.,62.,62.,62.,60.,62.]
maxft = [90.,90.,80.,80.,80.,77.,77.,80.,75.,80.,90.,82.,77.,80.]
minfgub = [70.,70.,60.,60.,70.,70.,60.,70.,70.,70.,70.,60.,60.,60.]
maxfgub = [70.,70.,80.,80.,70.,70.,90.,70.,70.,70.,82.,80.,80.,77.]
minfgmb = [70.,70.,60.,60.,70.,70.,62.,70.,70.,70.,70.,62.,62.,60.]
maxfgmb = [70.,70.,75.,75.,70.,70.,75.,70.,70.,70.,75.,80.,77.,77.]
minftb = [70.,70.,60.,60.,70.,70.,62.,70.,70.,70.,70.,62.,60.,62.]
maxftb = [70.,70.,77.,75.,70.,70.,75.,70.,70.,70.,75.,82.,80.,75.]
if full: 
    good_ind_t = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    good_ind_gu = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    good_ind_gm = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
if base_lim:
    good_ind_t = [0,1,2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20]
    good_ind_gu = [0,1,2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20]
    good_ind_gm = [0,1,2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20]
if gu_lim:
    good_ind_t = [0,1,2,3,4,5,6,8,9,10,12,15,17,18,19]
    good_ind_gu = [0,1,2,3,4,5,6,8,9,10,12,15,17,18,19]
    good_ind_gm = [0,1,2,3,4,5,6,8,9,10,12,15,17,18,19]
if gm_lim:
    good_ind_t = [0,2,3,4,5,8,9,10,13,14,15,16,17,18,19,20]
    good_ind_gu = [0,2,3,4,5,8,9,10,13,14,15,16,17,18,19,20]
    good_ind_gm = [0,2,3,4,5,8,9,10,13,14,15,16,17,18,19,20]
if all_type:
    good_ind_t = [0,1,2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20]
    good_ind_gu = [0,1,2,3,4,5,6,8,9,10,12,15,17,18,19]
    good_ind_gm = [0,2,3,4,5,8,9,10,13,14,15,16,17,18,19,20]

#Day 2 too small dataset and nans, Day 6 known bad data,
#Based on single day plots, setting updated datasets:
#if base_lim:
#    days = ['01','03','04','05','07','08','10','11','12','13','14']
#    two_art = [False,True,True,False,True,False,False,True,True,True,True]
#    count_a = [4196,3340,2506,1888,341,1738,1340,3182,2982,351,2684]
#    count_b = [0,3172,2325,0,133,0,0,1760,1504,233,1266]
#01 - good gmres, gures
#02 - too small
#03 - a, intermediate gm/gu, same with b (tres better than 1)
#04 - both pretty good
#05 - gm/t bad
#06 - bad noise data
#07 - a bad gu, okay gm, b good both
#08 - bad gm
#09 - all bad (extra noisy)
#10 - good gm/gu
#11 - a bad, b bad
#12 - a and b good
#13 - a okay, b bad
#14 - a okay, b bad
#GU optimized
#if gu_lim:
#    days = ['01','03','04','05','07','08','10','12','13','14']
#    two_art = [False,True,True,False,True,False,False,True,False,False]
#    count_a = [4196,3340,2506,1888,341,1738,1340,2982,351,2684]
#    count_b = [0,3172,2325,0,133,0,0,1504,0,0]
#GM optimized
#if gm_lim:
#    days = ['01','03','04','07','10','12','13','14']
#    two_art = [False,True,True,True,False,True,False,False]
#    count_a = [4196,3340,2506,341,1340,2982,351,2684]
#    count_b = [0,3172,2325,133,0,1504,0,0] 
bin = 9
tmean = []
gmmean = []
gumean = []
tres1 = []
tres2 = []
tres3 = []
tres4 = []
tres5 = []
tres6 = []
gmres1 = []
gmres2 = []
gmres3 = []
gmres4 = []
gmres5 = []
gmres6 = []
gures1 = []
gures2 = []
gures3 = []
gures4 = []
gures5 = []
gures6 = []
sub_names = []
sub_weights = []
tset = []
gmset = []
guset = []
subminft = []
submaxft = []
subminfgu = []
subminfgm = []
submaxfgu = []
submaxfgm = []
for d in range(0,len(days)):
    sub_names.append(days[d])
    sub_weights.append(count_a[d])
    if two_art[d]:
        single_name = days[d]+'b'
        sub_names.append(single_name)
        sub_weights.append(count_b[d])
    prefix = maindir+'June'+days[d]+'_day_to_night_'
    single_tmean = loadtxt(prefix+'tcal_mean.txt')
    single_gumean = loadtxt(prefix+'gucal_mean.txt') 
    single_gmmean = loadtxt(prefix+'gmcal_mean.txt') 
    if two_art[d]:
        single_tmeanb = loadtxt(prefix+'tcal_meanb.txt')
        single_gumeanb = loadtxt(prefix+'gucal_meanb.txt')
        single_gmmeanb = loadtxt(prefix+'gmcal_meanb.txt')

#    freqs = arange(40.6103764,128.8,1.2207528)
#    freqs = arange(40.122075279756,129.9,0.244150559512)
    freqs = arange(40.97660224,129.5,1.95320448)
#Setting up the data of interest for residual plot
    minfreq = where(freqs<=minft[int(days[d])-1])[0][-1]
    maxfreq = where(freqs<=maxft[int(days[d])-1])[0][-1]
    subminft.append(minfreq)
    submaxft.append(maxfreq)
    single_set = zeros(len(freqs))
    for i in range(minfreq,maxfreq):
        single_set[i] = 1.0
    tset.append(single_set)
    if two_art[d]:
        minfreqb = where(freqs<=minftb[int(days[d])-1])[0][-1]
        maxfreqb = where(freqs<=maxftb[int(days[d])-1])[0][-1]
        subminft.append(minfreqb)
        submaxft.append(maxfreqb)
        single_set = zeros(len(freqs))
        for i in range(minfreqb,maxfreqb):
            single_set[i] = 1.0 
        tset.append(single_set) 

    minfreq = where(freqs<=minfgu[int(days[d])-1])[0][-1]
    maxfreq = where(freqs<=maxfgu[int(days[d])-1])[0][-1]
    subminfgu.append(minfreq)
    submaxfgu.append(maxfreq)
    single_set = zeros(len(freqs))
    for i in range(minfreq,maxfreq):
        single_set[i] = 1.0 
    guset.append(single_set) 
    if two_art[d]:
        minfreqb = where(freqs<=minfgub[int(days[d])-1])[0][-1]
        maxfreqb = where(freqs<=maxfgub[int(days[d])-1])[0][-1]
        subminfgu.append(minfreqb)
        submaxfgu.append(maxfreqb) 
        single_set = zeros(len(freqs))
        for i in range(minfreqb,maxfreqb):
            single_set[i] = 1.0  
        guset.append(single_set)  

    minfreq = where(freqs<=minfgm[int(days[d])-1])[0][-1]
    maxfreq = where(freqs<=maxfgm[int(days[d])-1])[0][-1]
    subminfgm.append(minfreq)
    submaxfgm.append(maxfreq) 
    single_set = zeros(len(freqs))
    for i in range(minfreq,maxfreq):
        single_set[i] = 1.0 
    gmset.append(single_set) 
    if two_art[d]:
        minfreqb = where(freqs<=minfgmb[int(days[d])-1])[0][-1]
        maxfreqb = where(freqs<=maxfgmb[int(days[d])-1])[0][-1]
        single_set = zeros(len(freqs))
        subminfgm.append(minfreqb)
        submaxfgm.append(maxfreqb)  
        for i in range(minfreqb,maxfreqb):
            single_set[i] = 1.0  
        gmset.append(single_set)      

#    freqs = arange(41.09867752,128.,2.1973503)
#    print shape(freqs)
#    freqs = arange(40.122075279756,129.9,0.244150559512)
#    freqs = arange(40.1,130.1,90./len(single_tmean))
#    rb_tmean = []
#    rb_gmmean = []
#    rb_gumean = []
#    rb_freq = []
#    rb_tmeanb = []
#    rb_gumeanb = []
#    rb_gmmeanb = []
#    for i in range(0,len(single_tmean)/bin):
#        rb_tmean.append(ma.mean(single_tmean[bin*i:bin*(i+1)]))
#        rb_gmmean.append(ma.mean(single_gmmean[bin*i:bin*(i+1)]))
#        rb_gumean.append(ma.mean(single_gumean[bin*i:bin*(i+1)]))
#        rb_freq.append(ma.mean(freqs[bin*i:bin*(i+1)]))
#        if two_art[d]:
#            rb_tmeanb.append(ma.mean(single_tmeanb[bin*i:bin*(i+1)])) 
#            rb_gmmeanb.append(ma.mean(single_gmmeanb[bin*i:bin*(i+1)]))
#            rb_gumeanb.append(ma.mean(single_gumeanb[bin*i:bin*(i+1)]))
#
#    rb_tmean = array(rb_tmean)
#    rb_gmmean = array(rb_gmmean)
#    rb_gumean = array(rb_gumean)
#    if two_art[d]:
#        rb_tmeanb = array(rb_tmeanb)
#        rb_gmmeanb = array(rb_gmmeanb) 
#        rb_gumeanb = array(rb_gumeanb)
#    print shape(rb_tmeanb)

    rb_freq = freqs
    rb_tmean = single_tmean
    rb_gumean = single_gumean
    rb_gmmean = single_gmmean
    if two_art[d]:
        rb_tmeanb = single_tmeanb
        rb_gumeanb = single_gumeanb
        rb_gmmeanb = single_gmmeanb

    tmean.append(rb_tmean)
    gmmean.append(rb_gmmean)
    gumean.append(rb_gumean)
#    print shape(tmean)
    if two_art[d]:
        tmean.append(rb_tmeanb) 
        gmmean.append(rb_gmmeanb)
        gumean.append(rb_gumeanb)
#    print shape(tmean)

    for i in range(1,7):
        single_tres = loadtxt(prefix+'tcal_res'+str(i)+'.txt')
#        rb_single_tres = []
        single_gmres = loadtxt(prefix+'gmcal_res'+str(i)+'.txt')
#        rb_single_gmres = []
        single_gures = loadtxt(prefix+'gucal_res'+str(i)+'.txt')
#        rb_single_gures = []
        if two_art[d]:
            single_tresb = loadtxt(prefix+'tcal_res'+str(i)+'b.txt')
#            rb_single_tresb = []
            single_gmresb = loadtxt(prefix+'gmcal_res'+str(i)+'b.txt')
#            rb_single_gmresb = [] 
            single_guresb = loadtxt(prefix+'gucal_res'+str(i)+'b.txt')
#            rb_single_guresb = [] 
#        for j in range(0,len(single_tres)/bin):
#            print shape(single_tres[bin*j:bin*(j+1)])
#            rb_single_tres.append(ma.mean(single_tres[bin*j:bin*(j+1)]))
#            rb_single_gmres.append(ma.mean(single_gmres[bin*j:bin*(j+1)]))
#            rb_single_gures.append(ma.mean(single_gures[bin*j:bin*(j+1)]))
#            if two_art[d]:
#                rb_single_tresb.append(ma.mean(single_tresb[bin*j:bin*(j+1)]))
#                rb_single_gmresb.append(ma.mean(single_gmresb[bin*j:bin*(j+1)]))
#                rb_single_guresb.append(ma.mean(single_guresb[bin*j:bin*(j+1)]))             
#        rb_single_tres = array(rb_single_tres)
#        rb_single_gmres = array(rb_single_gmres)
#        rb_single_gures = array(rb_single_gures)
#        if two_art[d]:
#            rb_single_tresb = array(rb_single_tresb) 
#            rb_single_gmresb = array(rb_single_gmresb)
#            rb_single_guresb = array(rb_single_guresb)

        rb_single_tres = single_tres
        rb_single_gmres = single_gmres
        rb_single_gures = single_gures
        if two_art[d]:
            rb_single_tresb = single_tresb
            rb_single_gmresb = single_gmresb
            rb_single_guresb = single_guresb
        if i==1:
            tres1.append(rb_single_tres)
            gmres1.append(rb_single_gmres)
            gures1.append(rb_single_gures)
            if two_art[d]:
                tres1.append(rb_single_tresb)
                gmres1.append(rb_single_gmresb)
                gures1.append(rb_single_guresb)
        elif i==2:
            tres2.append(rb_single_tres)
            gmres2.append(rb_single_gmres)
            gures2.append(rb_single_gures)
            if two_art[d]:
                tres2.append(rb_single_tresb)
                gmres2.append(rb_single_gmresb)
                gures2.append(rb_single_guresb)
        elif i==3:
            tres3.append(rb_single_tres)
            gmres3.append(rb_single_gmres)
            gures3.append(rb_single_gures)
            if two_art[d]:
                tres3.append(rb_single_tresb)
                gmres3.append(rb_single_gmresb)
                gures3.append(rb_single_guresb)
        elif i==4:
            tres4.append(rb_single_tres)
            gmres4.append(rb_single_gmres)
            gures4.append(rb_single_gures)
            if two_art[d]:
                tres4.append(rb_single_tresb)
                gmres4.append(rb_single_gmresb)
                gures4.append(rb_single_guresb)
        elif i==5:
            tres5.append(rb_single_tres)
            gmres5.append(rb_single_gmres)
            gures5.append(rb_single_gures)
            if two_art[d]:
                tres5.append(rb_single_tresb)
                gmres5.append(rb_single_gmresb)
                gures5.append(rb_single_guresb)
        elif i==6:
            tres6.append(rb_single_tres)
            gmres6.append(rb_single_gmres)
            gures6.append(rb_single_gures)
            if two_art[d]:
                tres6.append(rb_single_tresb)
                gmres6.append(rb_single_gmresb)
                gures6.append(rb_single_guresb)

#print tres1[0][1:10]
print shape(tres1)
print shape(tmean)

rb_freq = array(rb_freq)
tmean = array(tmean)
gmmean = array(gmmean)
gumean = array(gumean)
tres1 = array(tres1)
tres2 = array(tres2)
tres3 = array(tres3)
tres4 = array(tres4)
tres5 = array(tres5)
tres6 = array(tres6)
gmres1 = array(gmres1)
gmres2 = array(gmres2)
gmres3 = array(gmres3)
gmres4 = array(gmres4)
gmres5 = array(gmres5)
gmres6 = array(gmres6)
gures1 = array(gures1)
gures2 = array(gures2)
gures3 = array(gures3)
gures4 = array(gures4)
gures5 = array(gures5)
gures6 = array(gures6)
tset = array(tset)
gmset = array(gmset)
guset = array(guset)
subminft = array(subminft)
submaxft = array(submaxft)
sub_weights = array(sub_weights)
print shape(subminft)


for i in range(0,len(sub_names)):
    pylab.plot(rb_freq,tmean[i],label='Resistor Cal')
    pylab.plot(rb_freq,gumean[i],label='GSM Cal')
    pylab.plot(rb_freq,gmmean[i],label='Mean Sub GSM Cal')
    pylab.grid()
    pylab.ylim(0,8e3)
    pylab.xlabel('Frequency (MHz)')
    pylab.xlim(60,90)
    pylab.legend()
    pylab.ylabel('Temperature (Kelvin)')
    pylab.title('Resistor Calibrated Daily Data')
    pylab.savefig(maindir+sub_names[i]+'_mean_comp',dpi=300)
    pylab.clf()
    

#for i in range(0,len(sub_names)):
#    pylab.plot(rb_freq,tmean[i],label=sub_names[i])
#pylab.grid()
#pylab.ylim(0,10e3)
#pylab.xlabel('Frequency (MHz)')
#pylab.xlim(60,90)
#pylab.legend()
#pylab.ylabel('Temperature (Kelvin)')
#pylab.title('Resistor Calibrated Daily Data')
#pylab.savefig(maindir+'tmean_comp',dpi=300)
#pylab.clf()

#for i in range(0,len(sub_names)):
#    pylab.plot(rb_freq,gumean[i],label=sub_names[i])
#pylab.grid()
#pylab.ylim(0,10e3)
#pylab.xlabel('Frequency (MHz)')
#pylab.xlim(60,90)
#pylab.ylabel('Temperature (Kelvin)')
#pylab.legend()
#pylab.title('GSM without Mean Sub Calibrated Daily Data')
#pylab.savefig(maindir+'gumean_comp',dpi=300)
#pylab.clf()

#for i in range(0,len(sub_names)):
#    pylab.plot(rb_freq,gmmean[i],label=sub_names[i])
#pylab.grid()
#pylab.ylim(0,10e3)
#pylab.xlabel('Frequency (MHz)')
#pylab.xlim(60,90)
#pylab.ylabel('Temperature (Kelvin)')
#pylab.legend()
#pylab.title('GSM with Mean Sub Calibrated Daily Data')
#pylab.savefig(maindir+'gmmean_comp',dpi=300)
#pylab.clf()

#for d in range(0,len(days)):
for d in range(0,len(sub_names)):
    pylab.plot(rb_freq,tres1[d],'b',label='n=1')
    pylab.plot(rb_freq,tres2[d],'g',label='n=2')
    pylab.plot(rb_freq,tres3[d],'c',label='n=3')
    pylab.plot(rb_freq,tres5[d],'r',label='n=5')
    pylab.axvspan(rb_freq[subminft[d]],rb_freq[submaxft[d]],facecolor='0.5',alpha=0.5)
#    pylab.vlines(rb_freq[subminft[d]],-50,50,'m')
#    pylab.vlines(rb_freq[submaxft[d]],-50,50,'m')
    pylab.ylim(-50,50)
    pylab.xlim(60,90)
    pylab.grid()
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Residuals (Kelvin)')
    pylab.legend()
    pylab.title('Resistor Calibrated Residuals')
    pylab.savefig(maindir+sub_names[d]+'_single_day_tres',dpi=300)
    pylab.clf()

for d in range(0,len(sub_names)):
    pylab.plot(rb_freq,gures1[d],label='n=1')
    pylab.plot(rb_freq,gures2[d],label='n=2')
    pylab.plot(rb_freq,gures3[d],label='n=3')
    pylab.plot(rb_freq,gures5[d],label='n=5')
    pylab.axvspan(rb_freq[subminfgu[d]],rb_freq[submaxfgu[d]],facecolor='0.5',alpha=0.5)
    pylab.ylim(-50,50)
    pylab.xlim(60,90)
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Residuals (Kelvin)')
    pylab.legend()
    pylab.grid()
    pylab.title('GSM without Mean Sub Calibrated Residuals')
    pylab.savefig(maindir+sub_names[d]+'_single_day_gures',dpi=300)
    pylab.clf()

for d in range(0,len(sub_names)):
    pylab.plot(rb_freq,gmres1[d],label='n=1')
    pylab.plot(rb_freq,gmres2[d],label='n=2')
    pylab.plot(rb_freq,gmres3[d],label='n=3')
    pylab.plot(rb_freq,gmres5[d],label='n=5')
    pylab.axvspan(rb_freq[subminfgm[d]],rb_freq[submaxfgm[d]],facecolor='0.5',alpha=0.5)
    pylab.ylim(-50,50)
    pylab.xlim(60,90)
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Residuals (Kelvin)')
    pylab.grid()
    pylab.legend()
    pylab.title('GSM with Mean Sub Calibrated Residuals')
    pylab.savefig(maindir+sub_names[d]+'_single_day_gmres',dpi=300)
    pylab.clf()

#print where(isnan(tres1))
#print shape(tres1)
#print shape(tmean)

# Separate out useable data:
tres2m = []
gures2m = []
gmres2m = []
tres2s = []
gures2s = []
gmres2s = []
for j in range(0,len(rb_freq)):
    single_tweight = []
    single_guweight = []
    single_gmweight = []
    single_tres = [] 
    single_gures = []
    single_gmres = []
    for g in range(0,len(good_ind_t)):
        i = good_ind_t[g]
        if tset[i,j] ==1.0:
            if abs(tres2[i,j])<50.: 
                single_tweight.append(sub_weights[i]) 
                single_tres.append(tres2[i,j])
    for g in range(0,len(good_ind_gu)):
        i = good_ind_gu[g]
        if guset[i,j] == 1.0:
            if abs(gures2[i,j])<50.: 
                single_guweight.append(sub_weights[i]) 
                single_gures.append(gures2[i,j])
    for g in range(0,len(good_ind_gm)):
        i = good_ind_gm[g]
        if gmset[i,j] == 1.0:
            if abs(gmres2[i,j])<50.: 
                single_gmweight.append(sub_weights[i]) 
                single_gmres.append(gmres2[i,j])
    if sum(tset[:,j])>1:
        single_tm = ma.average(single_tres,weights=single_tweight)
        single_ts = sqrt(ma.average((single_tres-single_tm)**2,weights=single_tweight))
        tres2m.append(single_tm)
        tres2s.append(single_ts)
    else: 
        tres2m.append(0.0)
        tres2s.append(0.0)
    if sum(guset[:,j])>1:
        single_gum = ma.average(single_gures,weights=single_guweight)
        single_gus = sqrt(ma.average((single_gures-single_gum)**2,weights=single_guweight))
        gures2m.append(single_gum)
        gures2s.append(single_gus)
    else:
        gures2m.append(0.0)
        gures2s.append(0.0) 
    if sum(gmset[:,j])>1:
        single_gmm = ma.average(single_gmres,weights=single_gmweight)
        single_gms = sqrt(ma.average((single_gmres-single_gmm)**2,weights=single_gmweight))
        gmres2m.append(single_gmm)
        gmres2s.append(single_gms)
    else: 
        gmres2m.append(0.0)
        gmres2s.append(0.0) 
    
#print tres2m

#print tres2s

#print rb_freq
#tres1m = ma.average(tres1,axis=0,weights=sub_weights)
#tres1s = sqrt(ma.average((tres1-tres1m)**2,axis=0,weights=sub_weights))
#tres2m = ma.average(tres2,axis=0,weights=sub_weights)
#tres2s = sqrt(ma.average((tres2-tres2m)**2,axis=0,weights=sub_weights))
#tres3m = ma.average(tres3,axis=0,weights=sub_weights)
#tres3s = sqrt(ma.average((tres3-tres3m)**2,axis=0,weights=sub_weights))
#tres4m = ma.average(tres4,axis=0,weights=sub_weights)
#tres4s = sqrt(ma.average((tres4-tres4m)**2,axis=0,weights=sub_weights))
#tres5m = ma.average(tres5,axis=0,weights=sub_weights)
#tres5s = sqrt(ma.average((tres5-tres5m)**2,axis=0,weights=sub_weights))
#tres6m = ma.average(tres6,axis=0,weights=sub_weights)
#tres6s = sqrt(ma.average((tres6-tres6m)**2,axis=0,weights=sub_weights))
#gures1m = ma.average(gures1,axis=0,weights=sub_weights)
#gures1s = sqrt(ma.average((gures1-gures1m)**2,axis=0,weights=sub_weights))
#gures2m = ma.average(gures2,axis=0,weights=sub_weights)
#gures2s = sqrt(ma.average((gures2-gures2m)**2,axis=0,weights=sub_weights))
#gures3m = ma.average(gures3,axis=0,weights=sub_weights)
#gures3s = sqrt(ma.average((gures3-gures3m)**2,axis=0,weights=sub_weights))
#gures4m = ma.average(gures4,axis=0,weights=sub_weights)
#gures4s = sqrt(ma.average((gures4-gures4m)**2,axis=0,weights=sub_weights))
#gures5m = ma.average(gures5,axis=0,weights=sub_weights)
#gures5s = sqrt(ma.average((gures5-gures5m)**2,axis=0,weights=sub_weights))
#gures6m = ma.average(gures6,axis=0,weights=sub_weights)
#gures6s = sqrt(ma.average((gures6-gures6m)**2,axis=0,weights=sub_weights))
#gmres1m = ma.average(gmres1,axis=0,weights=sub_weights)
#gmres1s = sqrt(ma.average((gmres1-gmres1m)**2,axis=0,weights=sub_weights))
#gmres2m = ma.average(gmres2,axis=0,weights=sub_weights)
#gmres2s = sqrt(ma.average((gmres2-gmres2m)**2,axis=0,weights=sub_weights))
#gmres3m = ma.average(gmres3,axis=0,weights=sub_weights)
#gmres3s = sqrt(ma.average((gmres3-gmres3m)**2,axis=0,weights=sub_weights))
#gmres4m = ma.average(gmres4,axis=0,weights=sub_weights)
#gmres4s = sqrt(ma.average((gmres4-gmres4m)**2,axis=0,weights=sub_weights))
#gmres5m = ma.average(gmres5,axis=0,weights=sub_weights)
#gmres5s = sqrt(ma.average((gmres5-gmres5m)**2,axis=0,weights=sub_weights))
#gmres6m = ma.average(gmres6,axis=0,weights=sub_weights)
#gmres6s = sqrt(ma.average((gmres6-gmres6m)**2,axis=0,weights=sub_weights))

#print tres1m
#print rb_freq
#First Plot how adding additional n's help for each model
#pylab.scatter(rb_freq,tres1m,c='b',edgecolor='b',s=5)
#pylab.errorbar(rb_freq,tres1m,tres1s,c='b',fmt=None)
pylab.scatter(rb_freq,tres2m,c='b',edgecolor='b',s=5)
pylab.errorbar(rb_freq,tres2m,tres2s,c='b',fmt=None)
#pylab.scatter(rb_freq,tres3m,c='g',edgecolor='g',s=5)
#pylab.errorbar(rb_freq,tres3m,tres3s,c='g',fmt=None)
#pylab.scatter(rb_freq,tres4m,c='c',edgecolor='c',s=5)
#pylab.errorbar(rb_freq,tres4m,tres4s,c='c',fmt=None)
#pylab.scatter(rb_freq,tres5m,c='m',edgecolor='m',s=5)
#pylab.errorbar(rb_freq,tres5m,tres5s,c='m',fmt=None)
#pylab.scatter(rb_freq,tres6m,c='k',edgecolor='k',s=5)
#pylab.errorbar(rb_freq,tres6m,tres6s,c='k',fmt=None)
pylab.xlim(60,90)
pylab.ylim(-50,50)
pylab.grid()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits: Res Cal Data')
if full:
    pylab.savefig(maindir+'tcal_res_mean_comp',dpi=300)
if base_lim:
    pylab.savefig(maindir+'tcal_res_mean_comp_baselim',dpi=300)
if gu_lim:
    pylab.savefig(maindir+'tcal_res_mean_comp_gulim',dpi=300)
if gm_lim:
    pylab.savefig(maindir+'tcal_res_mean_comp_gmlim',dpi=300)
pylab.clf()

#pylab.scatter(rb_freq,gures1m,c='b',edgecolor='b',s=5)
#pylab.errorbar(rb_freq,gures1m,gures1s,c='b',fmt=None)
pylab.scatter(rb_freq,gures2m,c='b',edgecolor='b',s=5)
pylab.errorbar(rb_freq,gures2m,gures2s,c='b',fmt=None)
#pylab.scatter(rb_freq,gures3m,c='r',edgecolor='r',s=5)
#pylab.errorbar(rb_freq,gures3m,gures3s,c='r',fmt=None)
#pylab.scatter(rb_freq,gures4m,c='c',edgecolor='c',s=5)
#pylab.errorbar(rb_freq,gures4m,gures4s,c='c',fmt=None)
#pylab.scatter(rb_freq,gures5m,c='m',edgecolor='m',s=5)
#pylab.errorbar(rb_freq,gures5m,gures5s,c='m',fmt=None)
#pylab.scatter(rb_freq,gures6m,c='k',edgecolor='k',s=5)
#pylab.errorbar(rb_freq,gures6m,gures6s,c='k',fmt=None)
pylab.xlim(60,90)
pylab.ylim(-10,10)
pylab.grid()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits: GSM Mean Unsub Data')
if full:
    pylab.savefig(maindir+'gucal_res_mean_comp',dpi=300)
if base_lim:
    pylab.savefig(maindir+'gucal_res_mean_comp_baselim',dpi=300)
if gu_lim:
    pylab.savefig(maindir+'gucal_res_mean_comp_gulim',dpi=300)
if gm_lim:
    pylab.savefig(maindir+'gucal_res_mean_comp_gmlim',dpi=300) 
#pylab.savefig(maindir+'gucal_res_mean_comp',dpi=300)
pylab.clf()

#pylab.scatter(rb_freq,gmres1m,c='b',edgecolor='b',s=5)
#pylab.errorbar(rb_freq,gmres1m,gmres1s,c='b',fmt=None)
pylab.scatter(rb_freq,gmres2m,c='b',edgecolor='b',s=5)
pylab.errorbar(rb_freq,gmres2m,gmres2s,c='b',fmt=None)
#pylab.scatter(rb_freq,gmres3m,c='r',edgecolor='r',s=5)
#pylab.errorbar(rb_freq,gmres3m,gmres3s,c='r',fmt=None)
#pylab.scatter(rb_freq,gmres4m,c='c',edgecolor='c',s=5)
#pylab.errorbar(rb_freq,gmres4m,gmres4s,c='c',fmt=None)
#pylab.scatter(rb_freq,gmres5m,c='m',edgecolor='m',s=5)
#pylab.errorbar(rb_freq,gmres5m,gmres5s,c='m',fmt=None)
#pylab.scatter(rb_freq,gmres6m,c='k',edgecolor='k',s=5)
#pylab.errorbar(rb_freq,gmres6m,gmres6s,c='k',fmt=None)
pylab.xlim(60,90) 
pylab.ylim(-20,20) 
pylab.grid() 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits: GSM Mean Subtract Data')
if full:
    pylab.savefig(maindir+'gmcal_res_mean_comp',dpi=300)
if base_lim: 
    pylab.savefig(maindir+'gmcal_res_mean_comp_baselim',dpi=300)
if gu_lim: 
    pylab.savefig(maindir+'gmcal_res_mean_comp_gulim',dpi=300)
if gm_lim: 
    pylab.savefig(maindir+'gmcal_res_mean_comp_gmlim',dpi=300)
#pylab.savefig(maindir+'gmcal_res_mean_comp',dpi=300)
pylab.clf()  

poslog_tres2m = []
poslog_tres2s1= []
poslog_tres2s2 = []
posf_tres2m = []
posarr_tres2m = []
posaf_tres2m = []
posdel_tres2m = []
neglog_tres2m = []
neglog_tres2s1= []
neglog_tres2s2 = []
negf_tres2m = []
negarr_tres2m = []
negaf_tres2m = []
negdel_tres2m = []
poslog_gures2m = []
poslog_gures2s1 = []
poslog_gures2s2 = []
posf_gures2m = []
posarr_gures2m = []
posaf_gures2m = []
posdel_gures2m = []
neglog_gures2m = []
neglog_gures2s1 = []
neglog_gures2s2 = []
negf_gures2m = []
negarr_gures2m = []
negaf_gures2m = []
negdel_gures2m = []
poslog_gmres2m = []
poslog_gmres2s1 = []
poslog_gmres2s2 = []
posf_gmres2m = []
posarr_gmres2m = []
posaf_gmres2m = []
posdel_gmres2m = []
neglog_gmres2m = []
neglog_gmres2s1 = []
neglog_gmres2s2 = []
negf_gmres2m = []
negarr_gmres2m = []
negaf_gmres2m = []
negdel_gmres2m = []

for i in range(0,len(rb_freq)):
    if tres2m[i]>0.:
        poslog_tres2m.append((tres2m[i]))
        if tres2m[i]-tres2s[i]>0.:
            poslog_tres2s1.append((tres2s[i]))
            poslog_tres2s2.append((tres2s[i]))
        else:
            poslog_tres2s1.append((0))
            poslog_tres2s2.append((tres2s[i]))
            posarr_tres2m.append(tres2m[i])
            posaf_tres2m.append(rb_freq[i])
            posdel_tres2m.append(0.5*tres2m[i])
        posf_tres2m.append(rb_freq[i])
    elif tres2m[i]<0.:
        neglog_tres2m.append((abs(tres2m[i])))
        if tres2m[i]+tres2s[i]<0.:
            neglog_tres2s1.append((tres2s[i]))
            neglog_tres2s2.append((tres2s[i]))
        else: 
            neglog_tres2s1.append((0))
            neglog_tres2s2.append((tres2s[i]))
            negarr_tres2m.append(abs(tres2m[i]))
            negaf_tres2m.append(rb_freq[i])
            negdel_tres2m.append(0.5*abs(tres2m[i]))
        negf_tres2m.append(rb_freq[i])
    if gures2m[i]>0.:
        poslog_gures2m.append((gures2m[i]))
        if gures2m[i]-gures2s[i]>0.:
            poslog_gures2s1.append((gures2s[i]))
            poslog_gures2s2.append((gures2s[i]))
        else:
            poslog_gures2s1.append((0))
            poslog_gures2s2.append((gures2s[i]))
            posarr_gures2m.append(gures2m[i])
            posaf_gures2m.append(rb_freq[i])
            posdel_gures2m.append(0.5*gures2m[i])
        posf_gures2m.append(rb_freq[i])
    elif gures2m[i]<0.:
        neglog_gures2m.append((abs(gures2m[i])))
        if gures2m[i]+gures2s[i]<0.:
            neglog_gures2s1.append((gures2s[i]))
            neglog_gures2s2.append((gures2s[i]))
        else:
            neglog_gures2s1.append((0))
            neglog_gures2s2.append((gures2s[i]))
            negarr_gures2m.append(abs(gures2m[i]))
            negaf_gures2m.append(rb_freq[i])
            negdel_gures2m.append(0.5*abs(gures2m[i]))
        negf_gures2m.append(rb_freq[i])
    if gmres2m[i]>0.:
        poslog_gmres2m.append((gmres2m[i]))
        if gmres2m[i]-gures2s[i]>0.:
            poslog_gmres2s1.append((gmres2s[i]))
            poslog_gmres2s2.append((gmres2s[i])) 
        else:
            poslog_gmres2s1.append((0))
            poslog_gmres2s2.append((gmres2s[i])) 
            posarr_gmres2m.append(gmres2m[i])
            posaf_gmres2m.append(rb_freq[i])
            posdel_gmres2m.append(0.5*gmres2m[i])
        posf_gmres2m.append(rb_freq[i])
    elif gmres2m[i]<0.:
        neglog_gmres2m.append((abs(gmres2m[i])))
        if gmres2m[i]+gmres2s[i]<0.:
            neglog_gmres2s1.append((gmres2s[i]))
            neglog_gmres2s2.append((gmres2s[i]))
        else:
            neglog_gmres2s1.append((0))
            neglog_gmres2s2.append((gmres2s[i])) 
            negarr_gmres2m.append(abs(gmres2m[i]))
            negaf_gmres2m.append(rb_freq[i])
            negdel_gmres2m.append(0.5*abs(gmres2m[i]))
        negf_gmres2m.append(rb_freq[i])

posf_tres2m = array(posf_tres2m)
negf_tres2m = array(negf_tres2m)
posf_gures2m = array(posf_gures2m)
negf_gures2m = array(negf_gures2m)
posf_gmres2m = array(posf_gmres2m)
negf_gmres2m = array(negf_gmres2m)
poslog_tres2s = array([poslog_tres2s1,poslog_tres2s2])
poslog_tres2m = array(poslog_tres2m)
neglog_tres2s = array([neglog_tres2s1,neglog_tres2s2])
neglog_tres2m = array(neglog_tres2m)
poslog_gures2s = array([poslog_gures2s1,neglog_gures2s2])
poslog_gures2m = array(poslog_gures2m)
neglog_gures2s = array([neglog_gures2s1,neglog_gures2s2])
neglog_gures2m = array(neglog_gures2m)
poslog_gmres2s = array([poslog_gmres2s1,poslog_gmres2s2])
poslog_gmres2m = array(poslog_gmres2m)
neglog_gmres2s = array([neglog_gmres2s1,neglog_gmres2s2])
neglog_gmres2m = array(neglog_gmres2m)
posarr_tres2m = array(posarr_tres2m)
negarr_tres2m = array(negarr_tres2m)
posaf_tres2m = array(posaf_tres2m)
negaf_tres2m = array(negaf_tres2m)
posdel_tres2m = array(posdel_tres2m)
negdel_tres2m = array(negdel_tres2m)
posdelf_tres2m = zeros(len(posaf_tres2m))
negdelf_tres2m = zeros(len(negaf_tres2m))
posarr_gures2m = array(posarr_gures2m)
negarr_gures2m = array(negarr_gures2m)
posaf_gures2m = array(posaf_gures2m)
negaf_gures2m = array(negaf_gures2m)
posdel_gures2m = array(posdel_gures2m)
negdel_gures2m = array(negdel_gures2m) 
posdelf_gures2m = zeros(len(posaf_gures2m))
negdelf_gures2m = zeros(len(negaf_gures2m))
posarr_gmres2m = array(posarr_gmres2m)
negarr_gmres2m = array(negarr_gmres2m)
posaf_gmres2m = array(posaf_gmres2m)
negaf_gmres2m = array(negaf_gmres2m)
posdel_gmres2m = array(posdel_gmres2m)
negdel_gmres2m = array(negdel_gmres2m) 
posdelf_gmres2m = zeros(len(posaf_gmres2m))
negdelf_gmres2m = zeros(len(negaf_gmres2m))
#print poslog_tres2m,neglog_tres2m
#print shape(poslog_tres2s),neglog_tres2s
afreq = anat_resid[:,0]
amean = anat_resid[:,1]
astd = anat_resid[:,2]
pos_af = []
neg_af = []
pos_am = []
neg_am = []
pos_as1 = []
pos_as2 = []
neg_as1 = []
neg_as2 = []
pos_arr = []
neg_arr = []
pos_arf = []
neg_arf = []
pos_del = []
neg_del = []
for i in range(0,len(afreq)):
    if amean[i]>0.:
        pos_af.append(afreq[i])
        pos_am.append(amean[i])
        if amean[i]-astd[i]>0.:
            pos_as1.append(astd[i])
            pos_as2.append(astd[i])
        else:
            pos_as1.append(0)
            pos_as2.append(astd[i])
            pos_arr.append(amean[i])
            pos_arf.append(afreq[i])
            pos_del.append(0.5*amean[i])
    else:
        neg_af.append(afreq[i])
        neg_am.append(abs(amean[i]))
        if amean[i]+astd[i]<0.:
            neg_as1.append(astd[i])
            neg_as2.append(astd[i])
        else:
            neg_as1.append(0)
            neg_as2.append(astd[i])
            neg_arr.append(abs(amean[i]))
            neg_arf.append(afreq[i])
            neg_del.append(0.5*abs(amean[i]))
pos_as = array([pos_as1,pos_as2])
neg_as = array([neg_as1,neg_as2])
pos_arr = array(pos_arr)
neg_arr = array(neg_arr)
pos_arf = array(pos_arf)
neg_arf = array(neg_arf)
pos_del = array(pos_del)
neg_del = array(neg_del)
pos_delf = zeros(len(pos_arf))
neg_delf = zeros(len(neg_arf))

th1f = th1[:,0]
th1m = abs(th1[:,3])
th2f = th2[:,0]
th2m = abs(th2[:,3])
th3f = th3[:,0]
th3m = abs(th3[:,3])

pylab.semilogy()
pylab.errorbar(posf_tres2m,poslog_tres2m,poslog_tres2s,ecolor='b',elw=8,fmt=None)
pylab.errorbar(negf_tres2m,neglog_tres2m,neglog_tres2s,ecolor='b',elw=8,fmt=None)
pylab.errorbar(posf_gmres2m-0.1,poslog_gmres2m,poslog_gmres2s,elw=8,ecolor='r',fmt=None)
pylab.errorbar(negf_gmres2m-0.1,neglog_gmres2m,neglog_gmres2s,elw=8,ecolor='r',fmt=None)
#pylab.errorbar(pos_af,pos_am,pos_as,elw=8,ecolor='g',fmt=None)
#pylab.errorbar(neg_af,neg_am,neg_as,elw=8,ecolor='g',fmt=None)
for i in range(0,len(posaf_tres2m)):
    pylab.annotate("",xy=(posaf_tres2m[i],posarr_tres2m[i]-posdel_tres2m[i]),xytext=(posaf_tres2m[i],posarr_tres2m[i]),arrowprops=dict(arrowstyle="->",color='b'),)
for i in range(0,len(negaf_tres2m)):
    pylab.annotate("",xy=(negaf_tres2m[i],negarr_tres2m[i]-negdel_tres2m[i]),xytext=(negaf_tres2m[i],negarr_tres2m[i]),arrowprops=dict(arrowstyle="->",color='b'),)
for i in range(0,len(posaf_gmres2m)):
    pylab.annotate("",xy=(posaf_gmres2m[i]-0.1,posarr_gmres2m[i]-posdel_gmres2m[i]),xytext=(posaf_gmres2m[i]-0.1,posarr_gmres2m[i]),arrowprops=dict(arrowstyle="->",color='r'),)
for i in range(0,len(negaf_gmres2m)):
    pylab.annotate("",xy=(negaf_gmres2m[i]-0.1,negarr_gmres2m[i]-negdel_gmres2m[i]),xytext=(negaf_gmres2m[i]-0.1,negarr_gmres2m[i]),arrowprops=dict(arrowstyle="->",color='r'),)
#for i in range(0,len(pos_arf)):
#    pylab.annotate("",xy=(pos_arf[i],pos_arr[i]-pos_del[i]),xytext=(pos_arf[i],pos_arr[i]),arrowprops=dict(arrowstyle="->",color='g'),)
#pylab.scatter(pos_af,pos_am,c='g',edgecolor='g',s=20,marker='o') 
#for i in range(0,len(neg_arf)):
#    pylab.annotate("",xy=(neg_arf[i],neg_arr[i]-neg_del[i]),xytext=(neg_arf[i],neg_arr[i]),arrowprops=dict(arrowstyle="->",color='g'),)
pylab.scatter(posf_tres2m,poslog_tres2m,c='b',edgecolor='b',s=30,marker='o')
pylab.scatter(negf_tres2m,neglog_tres2m,c='b',edgecolor='b',s=30,marker='s')
pylab.scatter(posf_gmres2m-0.1,poslog_gmres2m,c='r',edgecolor='r',s=30,marker='o')
pylab.scatter(negf_gmres2m-0.1,neglog_gmres2m,c='r',edgecolor='r',s=30,marker='s')
#pylab.scatter(pos_af,pos_am,c='g',edgecolor='g',s=30,marker='o') 
#pylab.scatter(neg_af,neg_am,c='g',edgecolor='g',s=30,marker='s') 
pylab.plot(th1f,th1m,'k-.')
pylab.plot(th2f,th2m,'k--')
pylab.plot(th3f,th3m,'k')
pylab.ylim(0.001,50)
pylab.xlim(60,90)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.grid()
if all_type:
    pylab.savefig(maindir+'joint_log_means',dpi=300)
pylab.clf()



