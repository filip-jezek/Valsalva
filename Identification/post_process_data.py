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


# CONTROL ROOM
DATA_FOLDER = R"..\data\Valsalva"
COST_FUNCTION_FOLDER = R"..\Identification\valsalva"
VALUE_LOG_DIRNAME = R'..\Schedules\\'
VALUE_LOG_FILENAME = '_current_costs.txt'
DRAW_PLOTS = True
# write the outputfiles
WRITE_FILE = False

# file_set = 'Sitting V_OO'
# files = ["V_00_sit_01", "V_00_sit_02", "V_00_sit_03", "V_00_sit_04"]
# file_set = 'All sitting'
# files = ["V_00_sit_01", "V_00_sit_02", "V_00_sit_03", "V_00_sit_04", "V_01_sit_01", "V_01_sit_02", "V_01_sit_03", "V_02_sit_01", "V_02_sit_02", "V_03_sit_01", "V_03_sit_02"]
# file_set = 'All supine'
# files = ["V_00_sup_01", "VEc_01_sup_01", "VEc_01_sup_02", "VEc_01_sup_03", "VEc_02_sup_01", "VEc_02_sup_02", "VEc_03_sup_01"]
file_set = 'mODEL'
files = 

# </ COTRNOL ROM

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

def writeTargetValues(targetValues, fileTag = ''):
    if not WRITE_FILE:
        return
    # print means and stds to file
    
    filename = 'targetValues' + fileTag + '.txt'
    
    with open(DATA_FOLDER + '\\' + filename, 'w') as file:
        file.write('Name, value, std\n')
        for name, targetValue in targetValues.items():
            value = numpy.mean(targetValue)
            std = numpy.std(targetValue)
            print("%s : %.4f +/- %.3f" % (name,value, std ))
            file.write('%s, %.4f, %.3f\n' % (name,value, std))

    print('--- Written to file %s ' % (DATA_FOLDER + '\\' + filename) )




def plotTargetValues(targetValues):
    # plot merged timecourse
    valsalva_start = 20
    valsalva_end = 35
    signal_end = 55
    baseline_bp = numpy.mean(targetValues['baseline_bp'])
    recovery_bp = numpy.mean(targetValues['ph5_recovery'])*baseline_bp
    baseline_hr = numpy.mean(targetValues['baseline_hr'])
    recovery_hr = numpy.mean(targetValues['ph5_hr_recovery'])*baseline_hr

    plt.figure()

    plt.plot((0, valsalva_start), [baseline_bp]*2, 'k')
    plt.errorbar((valsalva_start/2), baseline_bp, yerr=numpy.std(targetValues['baseline_bp']), fmt = 'b')
    plt.plot((signal_end - 5, signal_end), [recovery_bp]*2, 'k')
    plt.errorbar((signal_end - 2.5), recovery_bp, yerr=numpy.std(targetValues['ph5_hr_recovery'])*baseline_bp, fmt = 'b')

    plt.plot((0, valsalva_start), [baseline_hr]*2, 'g')
    plt.errorbar((valsalva_start/2), baseline_hr, yerr=numpy.std(targetValues['baseline_hr']), fmt = 'c')
    plt.plot((signal_end - 5, signal_end), [recovery_hr]*2, 'g')
    plt.errorbar((signal_end - 2.5), recovery_hr, yerr=numpy.std(targetValues['ph5_hr_recovery'])*baseline_bp, fmt = 'c')

    def plotMetric(t_val, t_offset, val, baseline, color):
        val_mean = numpy.mean(val)*baseline
        val_std = numpy.std(val)*baseline
        t_mean = numpy.mean(t_val) + t_offset
        t_std = numpy.std(t_val) 
        plt.errorbar(t_mean, val_mean, yerr= val_std, xerr=t_std, fmt = color)
        plt.plot(t_mean, val_mean, '*' + color)

    def plotBPMetric(t_val, t_offset, val, color = 'b'):
        plotMetric(t_val, t_offset, val, baseline_bp, color)

    def plotHRMetric(t_val, t_offset, val, color = 'c'):
        plotMetric(t_val, t_offset, val, baseline_hr, color)    

    plotBPMetric(targetValues['t_ph1_peak'], valsalva_start, targetValues['ph1_peak'])
    plotBPMetric(targetValues['t_ph2_mean_min'], valsalva_start, targetValues['ph2_mean_min'])
    plotBPMetric(targetValues['t_ph2_max'], valsalva_end, targetValues['ph2_max'])
    plotBPMetric(targetValues['t_ph4_drop'], valsalva_end, targetValues['ph4_drop'])
    plotBPMetric(targetValues['t_ph4_ovrshoot'], valsalva_end, targetValues['ph4_ovrshoot'])

    plotHRMetric(targetValues['t_ph1_hr_min'], valsalva_start, targetValues['ph1_hr_min'])
    plotHRMetric(targetValues['t_ph4_hr_max'], valsalva_end, targetValues['ph4_hr_max'])
    plotHRMetric(targetValues['t_ph4_hr_drop'], valsalva_end, targetValues['ph4_hr_drop'])
    plotHRMetric(targetValues['t_ph4_ovrshoot'], valsalva_end, targetValues['ph5_hr_recovery'])

    plt.title('Avg metrics of \'%s\' with %d elements' % (file_set, len(files)))

    plt.show(block = True)

targetValues = dict()
cf = importCostFunction()

for file in files:
    # filename = 'V_00_sit_01.mat'
    filename = file + '.mat'

    file_path = DATA_FOLDER + '\\' + filename
    matdata = skipy.loadmat(file_path)

    keyMapping = [ ['arterial_pressure', 'brachial_pressure'],
                ['heart_rate','heartRate.HR'],
                ['thoracic_pressure', 'thoracic_pressure']
            ]
    var_set = extractVars(matdata, keyMapping)

    if DRAW_PLOTS:
        var_set['__draw_plots'] = True
        var_set['__file_name'] = file
    
    var_set['__targetValuesFilename'] = 'targetValues_All_sitting.txt'

    objectives = cf.getObjectives(var_set)

    #print cost
    cost = sum(o.cost() for o in objectives)
    print("Cost of %s is %.4f" % (file, round(cost, 4)))


    for objective in objectives:
        if objective.name not in targetValues:
            targetValues[objective.name] = list()

        targetValues[objective.name].append(objective.value)

    # writeCost(objectives)

fileTag = '_' + file_set.replace(' ', '_')
writeTargetValues(targetValues, fileTag=fileTag)

plotTargetValues(targetValues)

pass
