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

VALUE_LOG_DIRNAME = '..\\Schedules\\'
VALUE_LOG_FILENAME =  '_current_costs.txt'
# // write the outputfiles
def writeCost(costs):
    total_cost = sum(costs)
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))

def logLine(val, cost, sum):
	# return ','.join(["%"val), str(cost), str(round(cost/sum*100))])
	return '%.3e,%.3e,%02d' % (val, cost, round(cost/sum*100))


def logOutput(vals, costs):
    # log the output, if the log directory exists. exit otherwise

	log_dirname = VALUE_LOG_DIRNAME
	# log_dirname = 'Schedules'

	if not os.path.isdir(log_dirname):
		log_dirname = '..\\'
		# cur_dirname = os.path.basename(os.getcwd())
		# log_filename = log_dirname + '\\' + cur_dirname + '_costs.txt'
	log_filename = log_dirname + VALUE_LOG_FILENAME
	with open(log_filename, 'a') as file:
		# prepare the line with value, cost value for this and percentage of total costs
		string_seq = map(lambda val, cost: logLine(val, cost, sum(costs)), vals, costs)
		file.write(', '.join(string_seq) + '\n')
	
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

vals = cf.getCurrentValues(var_set)

# check if file exists and build header otherwise
# no point in building header when we do not know the value names tho
if hasattr(cf, "getValueNames") and callable(cf.getValueNames):
	if not os.path.isfile(VALUE_LOG_DIRNAME + VALUE_LOG_FILENAME):
		with open(VALUE_LOG_DIRNAME + VALUE_LOG_FILENAME, 'w') as file:
			header_vals = cf.getValueNames()
			header = map(lambda s: s.rjust(9) + ',' + s.rjust(4) + '_cost, %', header_vals)
			file.write(', '.join(header) + '\n')
		# logOutput(cf.getCurrentValues(var_set), costs)

if hasattr(cf, "getCurrentValues") and callable(cf.getCurrentValues):
	logOutput(vals, costs)