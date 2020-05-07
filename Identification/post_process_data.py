import scipy.io as skipy
import numpy
import DyMat
# its imported dynamically based on current dir
# import cost_function as cf
import sys
import os
import time
import importlib.util
import re
from datetime import datetime
import fun_lib
import matplotlib.pyplot as plt
# from matplotlib import pyplot as plt
# import re
# import TerminalDS

DATA_FOLDER = R"..\data\Valsalva"
COST_FUNCTION_FOLDER = R"..\Identification\valsalva"
VALUE_LOG_DIRNAME = R'..\Schedules\\'
VALUE_LOG_FILENAME = '_current_costs.txt'
DRAW_PLOTS = True
# // write the outputfiles


def writeCost(objectives):

    total_cost = sum(o.cost() for o in objectives)
    # total_cost = sum(costs)
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))

def writeLogHeader(objectives):
    """check if file exists and build header otherwise
    """

    if not os.path.isfile(VALUE_LOG_DIRNAME + '\\' + VALUE_LOG_FILENAME):
        with open(VALUE_LOG_DIRNAME + '\\' + VALUE_LOG_FILENAME, 'w') as file:
            header = map(lambda o: o.name.rjust(5) + '_val,' + o.name.rjust(5) + '_trg, %', objectives)
            line = ',  '.join(header) + "  ,run, datetime"
            file.write(line + '\n')

def logLine(objective, total_cost):
    # return ','.join(["%"val), str(cost), str(round(cost/sum*100))])
    if objective.targetValue  is not None:
        return '%.3e,%.3e,%02d' % (objective.value, objective.targetValue , round(objective.cost()/total_cost*100))
    else:
        s = ' in limit' if objective.inLimit() else 'out limit'
        return '%.3e,%s,%02d' % (objective.value, s , round(objective.cost()/total_cost*100))


def logOutput(objectives):
    # log the output, if the log directory exists. exit otherwise

    log_dirname = VALUE_LOG_DIRNAME + '\\'
    run = 0
    # log_dirname = 'Schedules'

    if not os.path.isdir(log_dirname):
        log_dirname = '..\\'
    else:
        log_dirname = VALUE_LOG_DIRNAME + '\\'

    
    cur_dirname = os.path.basename(os.getcwd())
    run_match = re.match(r'[\w-]*-(\d+)$', cur_dirname)

    if run_match is not None:
        run = int(run_match[1])
    else:
        run = 0

        # log_filename = log_dirname + '\\' + cur_dirname + '_costs.txt'
    log_filename = log_dirname + VALUE_LOG_FILENAME
    with open(log_filename, 'a') as file:
        # prepare the line with value, cost value for this and percentage of total costs
        total_cost = sum(o.cost() for o in objectives)
        t = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        tail = "  ,%03d,%s" % (run, t)
        
        string_seq = map(lambda o: logLine(o, total_cost), objectives)

        file.write(',  '.join(string_seq) + tail + '\n')


def extractVars(matdata, keyMapping):

    matkeys, _ = zip(*keyMapping)
    
    var_set = {}
    # git first for the time
    time, _ = zip(*matdata[matkeys[0]])
    var_set['time'] = numpy.array(time)

    for matkey, datakey in keyMapping:
        _, data = zip(*matdata[matkey])
        var_set[datakey] = numpy.array(data)

    return var_set


def importCostFunction():

    spec = importlib.util.spec_from_file_location(
        'cost_function', COST_FUNCTION_FOLDER + '\\' + 'cost_function.py')
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)
    return cf


targetValues = dict()
cf = importCostFunction()

files = ["V_00_sit_01", "V_00_sit_02", "V_00_sit_03", "V_00_sit_04"]
for file in files:
    # filename = 'V_00_sit_01.mat'
    filename = file + '.mat'

    file_path = DATA_FOLDER + '\\' + filename
    matdata = skipy.loadmat(file_path)

    keyMapping = [ ['arterial_pressure', 'brachial_pressure'],
                ['heart_rate','heart_rate'],
                ['thoracic_pressure', 'thoracic_pressure']
            ]
    var_set = extractVars(matdata, keyMapping)

    if DRAW_PLOTS:
        var_set['__draw_plots'] = True
        plt.figure()
        plt.title(file)
        plt.show(block = False)

    objectives = cf.getObjectives(var_set)


    for objective in objectives:
        if objective.name not in targetValues:
            targetValues[objective.name] = list()

        targetValues[objective.name].append(objective.value)

    # writeCost(objectives)

for name, targetValue in targetValues.items():
    target_np = numpy.array(targetValue)
    print("%s : %.4f +/- %.3f" % (name,target_np.mean()*100, target_np.std()*100 ))

plt.show(block = True)
# writeLogHeader(objectives)
# logOutput(objectives)
