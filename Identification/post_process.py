import scipy.io as skipy
import DyMat
# its imported dynamically based on current dir
# import cost_function as cf
import sys
import os
import time
import importlib.util
# from matplotlib import pyplot as plt
# import re
# import TerminalDS

# // write the outputfiles
def writeCost(costs):
    total_cost = sum(costs)
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))

def logOutput(costs):
    # log the output, if the log directory exists. exit otherwise

	log_dirname = '..\\Schedules'
	# log_dirname = 'Schedules'

	if os.path.isdir(log_dirname):
		cur_dirname = os.path.basename(os.getcwd())
		log_filename = log_dirname + '\\' + cur_dirname + '_costs.txt'
		with open(log_filename, 'w') as file:
			file.write(str(costs))
	
def extractVars(d):
	# use just a subset or use all
	# vrs_names = ['Systemic1.posterior_tibial_T4_R236.u_C', 'Systemic1.aortic_arch_C2.port_a.pressure']
	# vrs_names = d.names(block = 2)
	vrs_names = d.names()
	timevar = d.abscissa(2)[0]

	var_set = {}
	var_set['time'] = timevar
	for vn in vrs_names:
		var_set[vn] = d.data(vn)
	return var_set		

def importCostFunction():

	spec = importlib.util.spec_from_file_location('cost_function', '..\\cost_function.py')
	cf = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(cf)
	return cf

tic = time.time()

filename = 'dsres.mat'
if len(sys.argv) > 1:
	filename = sys.argv[1]

d = DyMat.DyMatFile(filename)
toc = time.time()
print("Opening result in ", toc - tic, " s")

var_set = extractVars(d)

toc = time.time()
print("Loading result in ", toc - tic, " s")

cf = importCostFunction()
costs = cf.calculateCosts(var_set)

print("Calculating costs in ", time.time() - tic, " s")

writeCost(costs)
logOutput(costs)
