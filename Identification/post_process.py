import scipy.io as skipy
import DyMat
# its imported dynamically based on current dir
# import cost_function as cf
import sys
import os
import time
import importlib.util
import re
import fun_lib
from datetime import datetime
# from matplotlib import pyplot as plt
# import re
# import TerminalDS

VALUE_LOG_DIRNAME = '..\\Schedules'
VALUE_LOG_FILENAME = '_current_costs.txt'
DRAW_PLOTS = True
READ_OBJECTIVES = True
# If true loads ONLY non-constants
OMIT_LOADING_PARAMS = True
# // write the outputfiles


def writeCost(objectives):

    total_cost = fun_lib.countTotalSumCost(objectives)
    # total_cost = sum(costs)
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))
    print('Total costs: %s' % (total_cost))

def getLogFilePath():
    return fun_lib.getSafeLogDir(VALUE_LOG_DIRNAME) + VALUE_LOG_FILENAME

def writeLogHeader(objectives):
    """check if file exists and build header otherwise
    """
    filepath = getLogFilePath()
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as file:
            header = map(lambda o: o.name.rjust(5) + '_val,' + o.name.rjust(5) + '_trg, %', objectives)
            line = ',  '.join(header) + "  ,run, datetime, total_costs"
            file.write(line + '\n')

def logLine(objective : fun_lib.ObjectiveVar, total_cost):
    # return ','.join(["%"val), str(cost), str(round(cost/sum*100))])
    if objective.costFunctionType is fun_lib.CostFunctionType.DistanceFromZero:
        target = "%d" % 0
    elif objective.targetValue  is not None:
        target = "%.3e" % objective.targetValue
    else:
        target = ' in limit' if objective.inLimit() else 'out limit'
        
    return '%.3e,%s,%02d' % (objective.value, target , round(objective.cost()/total_cost*100))


def logOutput(objectives):
    # log the output, if the log directory exists. exit otherwise

    filepath = getLogFilePath()
    run = fun_lib.getRunNumber()
    with open(filepath, 'a') as file:
        # prepare the line with value, cost value for this and percentage of total costs
        total_cost = fun_lib.countTotalSumCost(objectives)
        t = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        tail = "  ,%03d,%s, %.6e" % (run, t, total_cost)
        
        string_seq = map(lambda o: logLine(o, total_cost), objectives)

        file.write(',  '.join(string_seq) + tail + '\n')


def extractVars(d):
    # use just a subset or use all
    # vrs_names = ['Systemic1.posterior_tibial_T4_R236.u_C', 'Systemic1.aortic_arch_C2.port_a.pressure']
    if OMIT_LOADING_PARAMS:
        vrs_names = d.names(block = 2)
    else:
        vrs_names = d.names()

    timevar = d.abscissa(2)[0]

    var_set = {}
    var_set['time'] = timevar
    for vn in vrs_names:
        var_set[vn] = d.data(vn)
    return var_set


def importCostFunction():

    spec = importlib.util.spec_from_file_location(
        'cost_function', '..\\cost_function.py')
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

if DRAW_PLOTS:
    var_set['__draw_plots'] = True
    var_set['__plot_title'] = "Run %i" % (fun_lib.getRunNumber())
    var_set['__saveFig_path'] = "%sFitFig_%03d.png" % (fun_lib.getSafeLogDir(VALUE_LOG_DIRNAME), fun_lib.getRunNumber())
    
objectives = cf.getObjectives(var_set)

print("Calculating costs in ", time.time() - tic, " s")

writeCost(objectives)

writeLogHeader(objectives)
logOutput(objectives)
