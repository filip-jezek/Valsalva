# reads the measured data and passes them to cost_function to calculate proper metrics. Metrics are averaged and saved to file for further use

import scipy.io as skipy
import numpy
import DyMat
# its imported dynamically based on selected dir
# import cost_function as cf
import sys
import os
import time
import importlib.util
import re
from datetime import datetime
import fun_lib
import matplotlib.pyplot as plt

from fun_lib import ObjectiveVar
# from matplotlib import pyplot as plt
# import re
# import TerminalDS


# CONTROL ROOM
DATA_FOLDER = R"..\data\Valsalva"
COST_FUNCTION_FOLDER = R"..\Identification\valsalva"
VALUE_LOG_DIRNAME = R'..\data\Valsalva\ProcessLog'
VALUE_LOG_FILENAME = '_current_costs.txt'

# draws plots for each individual file
DRAW_PLOTS = True
# write the outputfiles with targetValues
WRITE_FILE = True
READ_OBJECTIVES = True
USE_WEIGHING = True

# file_set = 'Sitting V_OO'
# files = ["V_00_sit_01", "V_00_sit_02", "V_00_sit_03", "V_00_sit_04"]
# file_set = 'All sitting'
# files = ["V_00_sit_01", "V_00_sit_02", "V_00_sit_03", "V_00_sit_04", "V_01_sit_01", "V_01_sit_02", "V_01_sit_03", "V_02_sit_01", "V_02_sit_02", "V_03_sit_01", "V_03_sit_02"]
file_set = 'All supine'
files = ["V_00_sup_01", "VEc_01_sup_01", "VEc_01_sup_02", "VEc_01_sup_03", "VEc_02_sup_01", "VEc_02_sup_02", "VEc_03_sup_01"]

# </ COTRNOL ROM

def processObjectiveMetrics(objectiveMetrics):
    """ Expects objectiveMetrics[subject][objective] = value and returns mean and std dicts
    """
    # get mean of the measurements for each subject
    subjectMeanMetrics = dict()
    for subject, subjectMetric in objectiveMetrics.items():
        for objective, objList in subjectMetric.items():
            # get mean of all the measurements per one subject
            if objective not in subjectMeanMetrics:
                subjectMeanMetrics[objective] = list()
            subjectMeanMetrics[objective].append(numpy.mean(objList))
    
    # get mean of all subject's measurements
    objectiveMeans = dict()
    for objective, objList in subjectMeanMetrics.items():
        objectiveMeans[objective] = numpy.mean(objList)

    # get variance per each subject
    subjectVariances = dict()
    for subject, subjectMetric in objectiveMetrics.items():
        for objective, objList in subjectMetric.items():
            objMean = objectiveMeans[objective]
            n_measurement = len(objList)
            if objective not in subjectVariances:
                subjectVariances[objective] = dict()    
            subjectVariances[objective][subject] = sum(1/n_measurement*(abs(x) - abs(objMean))**2 for x in objList)
    
    # sum the variances from all subjects
    objectiveStds = dict()
    for objective, subjectVariance in subjectVariances.items():
        n_subjects = len(subjectVariance)
        objectiveStds[objective] = numpy.sqrt(sum(x for x in subjectVariance.values())/n_subjects)
    
    # # prepare ttuples for the output
    # # std = sqrt(var)
    # targetValues = dict()
    # for key, objectiveMean in objectiveMeans.items():
    #     std = sqrt(objectiveVariance[key])
    #     targetValues[key] = (objectiveMean, std)

    return (objectiveMeans, objectiveStds)

        
def processCostMetrics(costMetrics : dict):
    for subject, cost_list in costMetrics.items():
        vals = (subject, numpy.mean(cost_list), numpy.std(cost_list))
        print('Subject %2d has mean cost %.4f with SD %.4f' % vals)


def getSubjectNum(s):
    # expecting the input as e.g. 'V_00_sup_01' it gets the second group, supposedly a participant number
    parts = s.split('_')
    return int(parts[1])

def getTargetFileName(file_set):
    fileTag = '_' + file_set.replace(' ', '_')
    return DATA_FOLDER + '\\' + 'targetValues' + fileTag + '.txt'

def getWeights(files):
    """ Finds occurences for the same filename root "AA_XX_bbb_" from AA_XX_bbb_YY to weight participants with multiple measurements
    """
    if not USE_WEIGHING:
        return [1]*len(files)

    weights = list()
    for file in files:
        count = sum(1 for x in files if file[0:-3] == x[0:-3])
        weights.append(1/count)
    return weights


def writeCost(objectives, filename='dsout.out'):

    total_cost = sum(o.cost() for o in objectives)
    # total_cost = sum(costs)
    with open(filename, 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))

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
    writeLogHeader(objectives)

    filepath = getLogFilePath()
    run = fun_lib.getRunNumber()
    with open(filepath, 'a') as file:
        # prepare the line with value, cost value for this and percentage of total costs
        total_cost = fun_lib.countTotalSumCost(objectives)
        t = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        tail = "  ,%03d,%s, %.6e" % (run, t, total_cost)
        
        wc = fun_lib.countTotalWeightedCost(objectives)
        string_seq = map(lambda o: logLine(o, wc), objectives)

        file.write(',  '.join(string_seq) + tail + '\n')

def extractVars(matdata, keyMapping):

    matkeys, _ = zip(*keyMapping)
    
    var_set = {}
    # git first for the time
    time, _ = zip(*matdata[matkeys[0]])
    var_set['time'] = numpy.array(time)

    for matkey, datakey in keyMapping:
        _, data = zip(*matdata[matkey])
        if data is not None:
            var_set[datakey] = numpy.array(data)

    return var_set


def importCostFunction():

    spec = importlib.util.spec_from_file_location(
        'cost_function', COST_FUNCTION_FOLDER + '\\' + 'cost_function.py')
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)
    return cf

def writeTargetValues(targetValues, targetStds, file_set):
    if not WRITE_FILE:
        return
    # print means and stds to file
    
    filename = getTargetFileName(file_set)
    with open(filename, 'w') as file:
        file.write('Name, value, std\n')
        for name in targetValues.keys():
            value = targetValues[name]
            std = targetStds[name]
            print("%s : %.4f +/- %.3f" % (name,value, std ))
            file.write('%s, %.4f, %.3f\n' % (name,value, std))

    print('--- Written to file %s ' % (filename) )




def plotTargetValues(targetValues, targetStds, styles = ('k', 'b', 'g', 'c')):
    """    PLots in the curent plot.
    styles = (bp_base_s, bp_s, hr_base_s, hr_s)
    """

    (bp_base_s, bp_s, hr_base_s, hr_s) = styles

    # plot merged timecourse
    valsalva_start = 20
    valsalva_end = 35
    signal_end = 55
    baseline_bp = targetValues['baseline_bp']
    recovery_bp = targetValues['ph5_recovery']*baseline_bp
    baseline_hr = targetValues['baseline_hr']
    recovery_hr = targetValues['ph5_hr_recovery']*baseline_hr

    plt.plot((0, valsalva_start), [baseline_bp]*2, bp_base_s)
    plt.errorbar((valsalva_start/2), baseline_bp, yerr=targetStds['baseline_bp'], fmt = bp_s)
    plt.plot((signal_end - 5, signal_end), [recovery_bp]*2, bp_base_s)
    plt.errorbar((signal_end - 2.5), recovery_bp, yerr=targetStds['ph5_hr_recovery']*baseline_bp, fmt = bp_s)

    plt.plot((0, valsalva_start), [baseline_hr]*2, hr_base_s)
    plt.errorbar((valsalva_start/2), baseline_hr, yerr=targetStds['baseline_hr'], fmt = hr_s)
    plt.plot((signal_end - 5, signal_end), [recovery_hr]*2, hr_base_s)
    plt.errorbar((signal_end - 2.5), recovery_hr, yerr=targetStds['ph5_hr_recovery']*baseline_bp, fmt = hr_s)

    def plotMetric(t_val, t_offset, val, baseline, color):
        val_mean = targetValues[val]*baseline
        val_std = targetStds[val]*baseline
        t_mean = targetValues[t_val] + t_offset
        t_std = targetStds[t_val]
        plt.errorbar(t_mean, val_mean, yerr= val_std, xerr=t_std, fmt = color)
        plt.plot(t_mean, val_mean, '*' + color)

    def plotBPMetric(t_val, t_offset, val, color = bp_s):
        plotMetric(t_val, t_offset, val, baseline_bp, color)

    def plotHRMetric(t_val, t_offset, val, color = hr_s):
        plotMetric(t_val, t_offset, val, baseline_hr, color)    

    plotBPMetric('t_ph1_peak', valsalva_start, 'ph1_peak')
    plotBPMetric('t_ph2_mean_min', valsalva_start, 'ph2_mean_min')
    plotBPMetric('t_ph2_max', valsalva_end, 'ph2_max')
    plotBPMetric('t_ph4_drop', valsalva_end, 'ph4_drop')
    plotBPMetric('t_ph4_ovrshoot', valsalva_end, 'ph4_ovrshoot')

    plotHRMetric('t_ph1_hr_min', valsalva_start, 'ph1_hr_min')
    plotHRMetric('t_ph4_hr_max', valsalva_end, 'ph4_hr_max')
    plotHRMetric('t_ph4_hr_drop', valsalva_end, 'ph4_hr_drop')
    plotHRMetric('t_ph4_ovrshoot', valsalva_end, 'ph5_hr_recovery')

    plt.title('Avg metrics of \'%s\' with %d elements' % (file_set, len(files)))
    plt.ylim(-10, 180)
    plt.xlim(0, 60)


objectiveMetrics = dict()
cf = importCostFunction()
# weights = getWeights(files)
costMetrics = dict()

for file in files:
    # filename = 'V_00_sit_01.mat'
    filename = file + '.mat'
    subject = getSubjectNum(file)

    file_path = DATA_FOLDER + '\\' + filename
    matdata = skipy.loadmat(file_path)

    keyMapping = [ ['arterial_pressure', 'brachial_pressure'],
                ['heart_rate','heartRate.HR'],
                ['thoracic_pressure', 'thoracic_pressure']
               # ['thoracic_pressure', 'Systemic1.thoracic_pressure']
            ]
    var_set = extractVars(matdata, keyMapping)

    if DRAW_PLOTS:
        var_set['__draw_plots'] = True

        var_set['__plot_title'] = file
        var_set['__saveFig_path'] = fun_lib.getSafeLogDir(COST_FUNCTION_FOLDER)  + file + '.png'
    
    if READ_OBJECTIVES:
       var_set['__targetValuesFilename'] = getTargetFileName(file_set)
    else:
        var_set['__targetValuesFilename'] = None
    

    objectives = cf.getObjectives(var_set)
    logOutput(objectives)

    for objective in objectives:
        if subject not in objectiveMetrics:
            objectiveMetrics[subject] = dict()
        if objective.name not in objectiveMetrics[subject]:
            objectiveMetrics[subject][objective.name] = list()

        objectiveMetrics[subject][objective.name].append(objective.value)
    
    if subject not in costMetrics:
        costMetrics[subject] = list()
    costMetrics[subject].append(fun_lib.countTotalWeightedCost(objectives))

    # writeCost(objectives)
(targetVals, targetStds) = processObjectiveMetrics(objectiveMetrics)
processCostMetrics(costMetrics)
writeTargetValues(targetVals, targetStds, file_set)

plt.figure()
plotTargetValues(targetVals, targetStds)
plt.savefig(getTargetFileName(file_set).replace('.txt', '.png') , dpi = 150)


pass

# compare sitting and supine targetvalues
def getAndPlotTargetVals(file_set, style):
    filename = getTargetFileName(file_set)
    objectives = fun_lib.updateObjectivesByValuesFromFile(filename)

    targetVals = dict()
    targetStds = dict()
    for o in objectives:
        targetVals[o.name] = o.targetValue
        targetStds[o.name] = o.std

    plotTargetValues(targetVals, targetStds, style)

plt.figure()
getAndPlotTargetVals('All sitting', ('k', 'b', 'g', 'c'))
getAndPlotTargetVals('All supine', ('k--', 'r', 'g--', 'm'))
plt.title('Sit (full, blue, cyan) vs. supine (dashed, red, magenta) comparisson')
plt.savefig(getTargetFileName('sit vs supine').replace('.txt', '.png') , dpi = 300)
