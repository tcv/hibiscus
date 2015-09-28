"""
Module to combine gsm data for all sidereal times into a single file for a given location and antenna.
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab
import os
import sys
sys.path.append(os.path.abspath('../../hibiscus'))

antenna = sys.argv[1]
indir = '../../supplemental_data/gsm_data_'+antenna+'_Karoo_test_update/'
outdir = indir
idate = '2015/4/1'
lon = '21.4109'
lat = '-30.7216'
elevation = 1080
gsm_freq = np.arange(50,111,1)

dirlist = os.listdir(indir)
full_data = np.zeros((len(dirlist),len(gsm_freq)))
full_times = np.zeros(len(dirlist))
tind= 0

for fname in dirlist:
    if fname.split('.')[-1]=='npy':
        if fname.split('_')[2]=='Karoo':
            time_label = fname.split('_')[-3]
            time = float(time_label.split('-')[0])+float(time_label.split('-')[1])/60.+float(time_label.split('-')[2])/3600.
            data = np.load(indir+fname)
            full_data[tind] = data
            full_times[tind] = time
            tind = tind+1

full_data = full_data[0:tind]
full_times = full_times[0:tind]
print min(full_times),max(full_times)

sortind = np.argsort(full_times)
sorttime = np.zeros(len(full_times))
sortdata = np.zeros((len(full_times),len(gsm_freq)))
for i in range(0,len(full_times)):
    sorttime[i] = full_times[sortind[i]]
    sortdata[i] = full_data[sortind[i]]

pylab.imshow(sortdata,vmin=0,vmax = 10000,aspect = 60./24., extent=(gsm_freq[0],gsm_freq[-1],sorttime[-1],sorttime[0]))
pylab.colorbar()
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Sidereal Time (Hours)')
pylab.title('Expected Temperature (Kelvin)')
pylab.savefig(outdir+'GSM_Temp_'+antenna+'_Karoo.png')
pylab.clf()

pylab.plot(gsm_freq,sortdata[0])
pylab.plot(gsm_freq,sortdata[len(sortdata)/4])
pylab.plot(gsm_freq,sortdata[len(sortdata)/2])
pylab.plot(gsm_freq,sortdata[3*len(sortdata)/4])
pylab.xlim(50,110)
pylab.ylim(0,6000)
pylab.xlabel('Frequency (MHz)')
pylab.ylabel('Temperature (Kelvin)')
pylab.savefig(outdir+'Sample_GSM_Temp_'+antenna+'_Karoo.png')
pylab.clf()


pylab.plot(sorttime,sortdata[:,20],label='70 MHz')
pylab.plot(sorttime,sortdata[:,30],label='80 MHz')
pylab.plot(sorttime,sortdata[:,40],label='90 MHz')
pylab.xlim(0,24.)
pylab.xlabel('Sidereal Time (Hours)')
pylab.ylim(0,5e3)
pylab.ylabel('Temperature (Kelvin)')
pylab.grid()
pylab.legend()
pylab.savefig(outdir+'Sample_GSM_Time_stream_'+antenna+'_Karoo.png')
pylab.clf()

np.save(outdir + 'gsm_data_full_'+antenna+'_Karoo.npy', sortdata)
np.save(outdir+'gsm_sid_time_full_'+antenna+'_Karoo.npy',sorttime)

