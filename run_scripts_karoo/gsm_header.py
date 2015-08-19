"""
Header for gsm generation code.
"""

import sys
import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy
import pylab
import scipy.interpolate as itp
import numpy.ma as ma
import scipy.optimize as opt
import os
import skrf as rf
import sys
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath('../../hibiscus'))
import file_funcs as ff
import gsm_funcs as gf
import gsm_data_generate as gsm
import gsm_generate_opt as gsm_opt
import time 

#This is the time in minutes being used for running the code.
myarg = sys.argv[1]
#myarg = arange(0,1440,3)
#myarg = [6,]
print myarg

#Other relevant information needed to run code
antenna = '100'
fplot = '150'

#Actual run
start = time.time()
gsm_opt.main(myarg,antenna,fplot)
finish=time.time()
print 'Code took ',finish-start,' seconds to run.'
