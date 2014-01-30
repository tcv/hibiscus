fitfunc = lambda p,x: (p[0]+p[1]*x)+sin(p[2]*2*pi*x+p[3])
errfunc = lambda p,x,y: fitfunc(p,x)-y
p0 = [-80,-0.02,0.01,-5]
#raw = short[-3]
#raw = raw[200:-400]
correct = []
for i in range(0,len(short_time)):
    raw = short[i]
    raw = raw[200:-400]
    result, success = optimize.leastsq(errfunc,p0[:],args=(new_freq[200:-400],raw),maxfev=1000)
    correct.append(sin(2*pi*new_freq*result[2]+result[3]))

diff = []
diff_time = []
short_corr = []
load_corr = []
for i in range(0,len(load_time)):
    for j in range(0,len(short_time)):
        if abs(load_time[i]-short_time[j])<0.005:
            diff.append(load[i]-short[j]+correct[j])
            diff_time.append(load_time[i])
            short_corr.append(short[j]-correct[j])
            load_corr.append(load[i])



correct_old =[]
for i in range(0,len(old_short_time)):
    raw = old_short[i]
    raw = raw[200:-400]
    result, success = optimize.leastsq(errfunc,p0[:],args=(new_freq[200:-400],raw),maxfev=1000)
    correct_old.append(sin(2*pi*new_freq*result[2]+result[3]))

diff_old = []
diff_time_old = []
short_corr_old = []
load_corr_old = []
for i in range(0,len(old_load_time)):
    for j in range(0,len(old_short_time)):
        if abs(old_load_time[i]-old_short_time[j])<0.005:
            diff_old.append(old_load[i]-old_short[j]+correct_old[j])
            diff_time_old.append(old_load_time[i])
            short_corr_old.append(old_short[j]-correct_old[j])
            load_corr_old.append(old_load[i])

old_gaincal_data = []
for i in range(0,len(old_data_effcal)):
    index = 0
    time_comp = 100.
    for j in range(0,len(diff_time_old)):
        if abs(diff_time_old[j]-old_ant_time[i])<time_comp:
            time_comp = abs(diff_time_old[j]-old_ant_time[i])
            index = j
    gain_cal_sm = fc.smooth(diff_old[index],new_freq,1.5)
    gain_corr_data = 300.*(old_data_effcal_db[i]-load_corr_old[index])/gain_cal_sm(new_freq)+300.
    old_gaincal_data.append(gain_corr_data)

median_old_gaincal = ma.median(old_gaincal_data,axis=0)

gaincal_data = []
for i in range(0,len(data_eff_cal)):
    index = 0
    time_comp = 100.
    for j in range(0,len(diff_time)):
        if abs(diff_time[j]-ant_time[i])<time_comp:
            time_comp = abs(diff_time[j]-ant_time[i])
            index = j
    gain_cal_sm = fc.smooth(diff[index],new_freq,1.5)
    gain_corr_data = 300.*(data_eff_cal_db[i]-load_corr[index])/gain_cal_sm(new_freq)+300.
    gaincal_data.append(gain_corr_data)

median_gaincal = ma.median(gaincal_data,axis=0)
