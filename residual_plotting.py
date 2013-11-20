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

maindir = '/lustre/tcv/mean_cal_data/'
#Day 2 too small dataset and nans, Day 6 known bad data,
#Based on single day plots, setting updated datasets:
#days = ['01','03','04','05','07','08','09','10','11','12','13','14']
#two_art = [False,True,True,False,True,False,False,False,True,True,True,True]
#count_a = [4196,3340,2506,1888,341,1738,3105,1340,3182,2982,351,2684]
#count_b = [0,3172,2325,0,133,0,0,0,1760,1504,233,1266]
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
#days = ['01','03','04','05','07','08','10','12','13','14']
#two_art = [False,True,True,False,True,False,False,True,False,False]
#count_a = [4196,3340,2506,1888,341,1738,1340,2982,351,2684]
#count_b = [0,3172,2325,0,133,0,0,1504,0,0]
#GM optimized
days = ['01','03','04','07','10','12','13','14']
two_art = [False,True,True,True,False,True,False,False]
count_a = [4196,3340,2506,341,1340,2982,351,2684]
count_b = [0,3172,2325,133,0,1504,0,0] 
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
    freqs = arange(40.122075279756,129.9,0.244150559512)
#    freqs = arange(40.1,130.1,90./len(single_tmean))
    rb_tmean = []
    rb_gmmean = []
    rb_gumean = []
    rb_freq = []
    rb_tmeanb = []
    rb_gumeanb = []
    rb_gmmeanb = []
    for i in range(0,len(single_tmean)/bin):
        rb_tmean.append(ma.mean(single_tmean[bin*i:bin*(i+1)]))
        rb_gmmean.append(ma.mean(single_gmmean[bin*i:bin*(i+1)]))
        rb_gumean.append(ma.mean(single_gumean[bin*i:bin*(i+1)]))
        rb_freq.append(ma.mean(freqs[bin*i:bin*(i+1)]))
        if two_art[d]:
            rb_tmeanb.append(ma.mean(single_tmeanb[bin*i:bin*(i+1)])) 
            rb_gmmeanb.append(ma.mean(single_gmmeanb[bin*i:bin*(i+1)]))
            rb_gumeanb.append(ma.mean(single_gumeanb[bin*i:bin*(i+1)]))

    rb_tmean = array(rb_tmean)
    rb_gmmean = array(rb_gmmean)
    rb_gumean = array(rb_gumean)
    if two_art[d]:
        rb_tmeanb = array(rb_tmeanb)
        rb_gmmeanb = array(rb_gmmeanb) 
        rb_gumeanb = array(rb_gumeanb)
#    print shape(rb_tmeanb)

#    rb_freq = freqs
#    rb_tmean = single_tmean
#    rb_gumean = single_gumean
#    rb_gmmean = single_gmmean
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
        rb_single_tres = []
        single_gmres = loadtxt(prefix+'gmcal_res'+str(i)+'.txt')
        rb_single_gmres = []
        single_gures = loadtxt(prefix+'gucal_res'+str(i)+'.txt')
        rb_single_gures = []
        if two_art[d]:
            single_tresb = loadtxt(prefix+'tcal_res'+str(i)+'b.txt')
            rb_single_tresb = []
            single_gmresb = loadtxt(prefix+'gmcal_res'+str(i)+'b.txt')
            rb_single_gmresb = [] 
            single_guresb = loadtxt(prefix+'gucal_res'+str(i)+'b.txt')
            rb_single_guresb = [] 
        for j in range(0,len(single_tres)/bin):
#            print shape(single_tres[bin*j:bin*(j+1)])
            rb_single_tres.append(ma.mean(single_tres[bin*j:bin*(j+1)]))
            rb_single_gmres.append(ma.mean(single_gmres[bin*j:bin*(j+1)]))
            rb_single_gures.append(ma.mean(single_gures[bin*j:bin*(j+1)]))
            if two_art[d]:
                rb_single_tresb.append(ma.mean(single_tresb[bin*j:bin*(j+1)]))
                rb_single_gmresb.append(ma.mean(single_gmresb[bin*j:bin*(j+1)]))
                rb_single_guresb.append(ma.mean(single_guresb[bin*j:bin*(j+1)]))             
        rb_single_tres = array(rb_single_tres)
        rb_single_gmres = array(rb_single_gmres)
        rb_single_gures = array(rb_single_gures)
        if two_art[d]:
            rb_single_tresb = array(rb_single_tresb) 
            rb_single_gmresb = array(rb_single_gmresb)
            rb_single_guresb = array(rb_single_guresb)

#        rb_single_tres = single_tres
#        rb_single_gmres = single_gmres
#        rb_single_gures = single_gures
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

for i in range(0,len(sub_names)):
    pylab.plot(rb_freq,tmean[i],label=sub_names[i])
pylab.grid()
pylab.ylim(0,10e3)
pylab.xlabel('Frequency (MHz)')
pylab.xlim(60,90)
pylab.legend()
pylab.ylabel('Temperature (Kelvin)')
pylab.title('Resistor Calibrated Daily Data')
pylab.savefig(maindir+'tmean_comp',dpi=300)
pylab.clf()

for i in range(0,len(sub_names)):
    pylab.plot(rb_freq,gumean[i],label=sub_names[i])
pylab.grid()
pylab.ylim(0,10e3)
pylab.xlabel('Frequency (MHz)')
pylab.xlim(60,90)
pylab.ylabel('Temperature (Kelvin)')
pylab.legend()
pylab.title('GSM without Mean Sub Calibrated Daily Data')
pylab.savefig(maindir+'gumean_comp',dpi=300)
pylab.clf()

for i in range(0,len(sub_names)):
    pylab.plot(rb_freq,gmmean[i],label=sub_names[i])
pylab.grid()
pylab.ylim(0,10e3)
pylab.xlabel('Frequency (MHz)')
pylab.xlim(60,90)
pylab.ylabel('Temperature (Kelvin)')
pylab.legend()
pylab.title('GSM with Mean Sub Calibrated Daily Data')
pylab.savefig(maindir+'gmmean_comp',dpi=300)
pylab.clf()

#for d in range(0,len(days)):
for d in range(0,len(sub_names)):
    pylab.plot(rb_freq,tres1[d],label='n=1')
    pylab.plot(rb_freq,tres2[d],label='n=2')
    pylab.plot(rb_freq,tres3[d],label='n=3')
    pylab.plot(rb_freq,tres6[d],label='n=6')
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
    pylab.plot(rb_freq,gures6[d],label='n=6')
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
    pylab.plot(rb_freq,gmres6[d],label='n=6')
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


tres1m = ma.average(tres1,axis=0,weights=sub_weights)
tres1s = sqrt(ma.average((tres1-tres1m)**2,axis=0,weights=sub_weights))
tres2m = ma.average(tres2,axis=0,weights=sub_weights)
tres2s = sqrt(ma.average((tres2-tres2m)**2,axis=0,weights=sub_weights))
tres3m = ma.average(tres3,axis=0,weights=sub_weights)
tres3s = sqrt(ma.average((tres3-tres3m)**2,axis=0,weights=sub_weights))
tres4m = ma.average(tres4,axis=0,weights=sub_weights)
tres4s = sqrt(ma.average((tres4-tres4m)**2,axis=0,weights=sub_weights))
tres5m = ma.average(tres5,axis=0,weights=sub_weights)
tres5s = sqrt(ma.average((tres5-tres5m)**2,axis=0,weights=sub_weights))
tres6m = ma.average(tres6,axis=0,weights=sub_weights)
tres6s = sqrt(ma.average((tres6-tres6m)**2,axis=0,weights=sub_weights))
gures1m = ma.average(gures1,axis=0,weights=sub_weights)
gures1s = sqrt(ma.average((gures1-gures1m)**2,axis=0,weights=sub_weights))
gures2m = ma.average(gures2,axis=0,weights=sub_weights)
gures2s = sqrt(ma.average((gures2-gures2m)**2,axis=0,weights=sub_weights))
gures3m = ma.average(gures3,axis=0,weights=sub_weights)
gures3s = sqrt(ma.average((gures3-gures3m)**2,axis=0,weights=sub_weights))
gures4m = ma.average(gures4,axis=0,weights=sub_weights)
gures4s = sqrt(ma.average((gures4-gures4m)**2,axis=0,weights=sub_weights))
gures5m = ma.average(gures5,axis=0,weights=sub_weights)
gures5s = sqrt(ma.average((gures5-gures5m)**2,axis=0,weights=sub_weights))
gures6m = ma.average(gures6,axis=0,weights=sub_weights)
gures6s = sqrt(ma.average((gures6-gures6m)**2,axis=0,weights=sub_weights))
gmres1m = ma.average(gmres1,axis=0,weights=sub_weights)
gmres1s = sqrt(ma.average((gmres1-gmres1m)**2,axis=0,weights=sub_weights))
gmres2m = ma.average(gmres2,axis=0,weights=sub_weights)
gmres2s = sqrt(ma.average((gmres2-gmres2m)**2,axis=0,weights=sub_weights))
gmres3m = ma.average(gmres3,axis=0,weights=sub_weights)
gmres3s = sqrt(ma.average((gmres3-gmres3m)**2,axis=0,weights=sub_weights))
gmres4m = ma.average(gmres4,axis=0,weights=sub_weights)
gmres4s = sqrt(ma.average((gmres4-gmres4m)**2,axis=0,weights=sub_weights))
gmres5m = ma.average(gmres5,axis=0,weights=sub_weights)
gmres5s = sqrt(ma.average((gmres5-gmres5m)**2,axis=0,weights=sub_weights))
gmres6m = ma.average(gmres6,axis=0,weights=sub_weights)
gmres6s = sqrt(ma.average((gmres6-gmres6m)**2,axis=0,weights=sub_weights))

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
pylab.ylim(-100,100)
pylab.grid()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits: Res Cal Data')
pylab.savefig(maindir+'tcal_res_mean_comp_gmopt',dpi=300)
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
pylab.ylim(-5,5)
pylab.grid()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits: GSM Mean Unsub Data')
pylab.savefig(maindir+'gucal_res_mean_comp_gmopt',dpi=300)
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
pylab.ylim(-50,50) 
pylab.grid() 
pylab.xlabel('Frequency (MHz)') 
pylab.ylabel('Residuals (Kelvin)')
pylab.title('Residuals from log fits: GSM Mean Subtract Data')
pylab.savefig(maindir+'gmcal_res_mean_comp_gmopt',dpi=300)
pylab.clf()  

#pylab.savefig('test.png')
