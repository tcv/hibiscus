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

#myarg=int(sys.argv[1])
#print "python thinks your arg is " + repr(myarg) + " which when doubled is " + repr(2*myarg)

myarg = sys.argv[1]
print myarg
gsm.main(myarg)
