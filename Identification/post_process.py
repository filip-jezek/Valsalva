import scipy.io as skipy
import DyMat
import cost_function as cf
import sys
import time
# from matplotlib import pyplot as plt
# import re
# import TerminalDS


tic = time.time()

filename = 'dsres.mat'
if len(sys.argv) > 1:
	filename = sys.argv[1]

d = DyMat.DyMatFile(filename)
timevar = d.abscissa(2)[0]

toc = time.time()
print("Loading result in ", toc - tic, " s")

# use just a subset or use all
# vrs_names = ['Systemic1.posterior_tibial_T4_R236.u_C', 'Systemic1.aortic_arch_C2.port_a.pressure']
# vrs_names = d.names(block = 2)
vrs_names = d.names()


var_set = {}
var_set['time'] = timevar
for vn in vrs_names:
	var_set[vn] = d.data(vn)

cf.WriteCost(var_set)
print("Calculating costs in ", time.time() - toc, " s")