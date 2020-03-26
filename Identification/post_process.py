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
VALUE_LOG_FILENAME = '_current_costs.txt'
# // write the outputfiles


def writeCost(objectives):

    total_cost = sum(o.cost() for o in objectives)
    # total_cost = sum(costs)
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))


def logLine(objective, total_cost):
    # return ','.join(["%"val), str(cost), str(round(cost/sum*100))])
    if objective.targetValue  is not None:
        return '%.3e,%.3e,%02d' % (objective.value, objective.targetValue , round(objective.cost()/total_cost*100))
    else:
        s = ' in limit' if objective.inLimit() else 'out limit'
        return '%.3e,%s,%02d' % (objective.value, s , round(objective.cost()/total_cost*100))


def logOutput(objectives):
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
        total_cost = sum(o.cost() for o in objectives)
        string_seq = map(lambda o: logLine(o, total_cost), objectives)

        file.write(',  '.join(string_seq) + '\n')


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
objectives = cf.getObjectives(var_set)

print("Calculating costs in ", time.time() - tic, " s")

writeCost(objectives)

# check if file exists and build header otherwise
# no point in building header when we do not know the value names tho
if not os.path.isfile(VALUE_LOG_DIRNAME + VALUE_LOG_FILENAME):
    with open(VALUE_LOG_DIRNAME + VALUE_LOG_FILENAME, 'w') as file:
        header = map(lambda o: o.name.rjust(5) + '_val,' + o.name.rjust(5) + '_trg, %', objectives)
        file.write(',  '.join(header) + '\n')
    # logOutput(cf.getCurrentValues(var_set), costs)

logOutput(objectives)
