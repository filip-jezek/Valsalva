import scipy.io as skipy
import DyMat
# its imported dynamically based on current dir
# import cost_function as cf
import sys
import os
import time
import re
import fun_lib
from datetime import datetime
import json
# import ast
# from matplotlib import pyplot as plt
# import re
# import TerminalDS

# SETTINGSFILE = 'post_process_options.json'
# // write the outputfiles


def writeCost(objectives):

    total_cost = fun_lib.countTotalSumCost(objectives)
    # total_cost = sum(costs)
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + repr(total_cost))
    print('Total costs: %s' % (total_cost))

def getLogFilePath(var_set):
    return fun_lib.getSafeLogDir(var_set['__VALUE_LOG_DIRNAME']) + var_set['__VALUE_LOG_FILENAME']

def writeLogHeader(objectives, var_set):
    """check if file exists and build header otherwise
    """
    filepath = getLogFilePath(var_set)
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as file:
            header = map(lambda o: o.name.rjust(5) + '_val,' + o.name.rjust(5) + '_trg, %', objectives)
            line = ',  '.join(header) + "  ,run, datetime, total_costs"
            file.write(line + '\n')

def logLine(objective : fun_lib.ObjectiveVar, total_cost):
    # return ','.join(["%"val), str(cost), str(round(cost/sum*100))])

        
    return '%.3e,%s,%02d' % (objective.value, objective.target_log() , round(objective.cost()/total_cost*100))


def logOutput(objectives, var_set):
    # log the output, if the log directory exists. exit otherwise

    filepath = getLogFilePath(var_set)
    run = fun_lib.getRunNumber()
    with open(filepath, 'a') as file:
        # prepare the line with value, cost value for this and percentage of total costs
        total_cost = fun_lib.countTotalSumCost(objectives)
        t = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        tail = "  ,%03d,%s, %.6e" % (run, t, total_cost)
        
        tc = fun_lib.countTotalSumCost(objectives)
        string_seq = map(lambda o: logLine(o, tc), objectives)

        file.write(',  '.join(string_seq) + tail + '\n')


def extractVars(d, var_set, prefix = None):
    # use just a subset or use all
    # vrs_names = ['Systemic1.posterior_tibial_T4_R236.u_C', 'Systemic1.aortic_arch_C2.port_a.pressure']
    if var_set['__OMIT_LOADING_PARAMS']:
        vrs_names = d.names(block = 2)
    else:
        vrs_names = d.names()

    timevar = d.abscissa(2)[0]

    if prefix is None or prefix == '__None':
        var_set['time'] = timevar
    else:
        var_set[prefix + '.' +  'time'] = timevar

    for vn in vrs_names:
        if prefix is None or prefix == '__None':
            var_set[vn] = d.data(vn)
        else:
            var_set[prefix + '.' + vn] = d.data(vn)
           
    return var_set


def loadInputFiles(var_set):
    
    def loadMatFile(filename, prefix = None):
        tic = time.time()
        d = DyMat.DyMatFile(filename)
        toc = time.time()
        print("Opening %s in %.3f ms" % (filename, (toc - tic)*1000))

        extractVars(d, var_set, prefix = prefix)

    if '__inputFiles' in var_set and var_set['__inputFiles'] is not None:
        inputFile_args = var_set['__inputFiles']
    else:
        inputFile_args = sys.argv[1:]

    if len(inputFile_args) == 0:
        return loadMatFile('dsres.mat')
    elif len(inputFile_args) == 1:
        return loadMatFile(inputFile_args[0])
    elif len(inputFile_args) > 2:
        # multiple files as parameters in the form:
        # file1 prefix1 file2 prefix2 file3 prefix3
        if len(inputFile_args) % 2 == 1:
            raise ValueError("The argument list must be even, i.e. file1 prefix1 file2 prefix2 file3 prefix3. Use __None for no prefix")
        files = inputFile_args[::2]
        prefixes = inputFile_args[1::2]

        for f, p in zip(files, prefixes):
            loadMatFile(f, p)
        
        # var_setUpdate = {k: v for d in var_set_gen for k, v in d.items()}
        # var_set.update(var_setUpdate)

        # get current dir as a case name
        var_set['__plot_title'] = "Case %s" % (os.path.basename(os.getcwd()))
        var_set['__saveFig_path'] = "FitFig_0.png"

        return var_set

def getObjectives(var_set) -> fun_lib.ObjectiveVar:
    tic = time.time()
    cf = fun_lib.importCostFunction(dir = var_set['__COST_FUNC_PATH'])

    if var_set['__DRAW_PLOTS_OVERRIDE']:
        var_set['__draw_plots'] = True
        var_set['__plot_title'] = "Run %i" % (fun_lib.getRunNumber())
        var_set['__saveFig_path'] = "%sFitFig_%03d.png" % (fun_lib.getSafeLogDir(var_set['__VALUE_LOG_DIRNAME']), fun_lib.getRunNumber())
        
    # log combined output if supported by the cost function
    var_set['__objectivesLog_path'] = "%sobjectiveLog_%03d.csv" % (fun_lib.getSafeLogDir(var_set['__VALUE_LOG_DIRNAME']), fun_lib.getRunNumber())

    objectives = cf.getObjectives(var_set)

    print("Calculating costs in ", time.time() - tic, " s")
    return objectives

def logCrash(logdirname, line:str):
    with open(fun_lib.getSafeLogDir(logdirname) + 'errorLog.txt', 'a') as f:
        s = '%s: At run %d sumfin went wong: %s\n' % (datetime.now(), fun_lib.getRunNumber(), line)
        print(s)
        f.write(s)
        
def importOptions(var_set:dict, settingsFile = 'post_process_options.json'):
    """ Updates params from settings file (returns true) or loads default (returns False). 
    """
    defaultOptionsString = r"""{ 
        "__VALUE_LOG_DIRNAME" : "..\\Schedules",
        "__VALUE_LOG_FILENAME" : "_current_costs.txt",
        "__DRAW_PLOTS_OVERRIDE" : true,
        "__READ_OBJECTIVES" : true,
        "__OMIT_LOADING_PARAMS" : true,
		"__COST_FUNC_PATH" : "..\\",
		"__SORT_COSTS_BY" : "costs",
		"__LOG_OBJECTIVES_PATH" : "..\\Schedules\\objectives_log%03d.txt",
        "__LOG_OBJECTIVES_PATH_BASE" : null,
		"__LOG_ALL_OBJECTIVES_PATH" : null,
        "__LOG_ALL_OBJECTIVES_PATH_BASE" : null,
		"__targetValuesFilename" : "../../../data/Valsalva/targetValues_All_supine.txt",
		"__inputFiles" : ["dsres.mat"],
		"__root_path" : "..\\..\\"

}"""

    if not os.path.exists(settingsFile):
        sets = json.loads(defaultOptionsString)
        var_set.update(sets) 
        return False   
    else:
        with open(settingsFile) as file:
            sets = json.load(file)
            var_set.update(sets)
            return True

def processDyMatFile():
    
    var_set = {}
    # import options
    importOptions(var_set)

    try:
        
        loadInputFiles(var_set)

        objectives = getObjectives(var_set)

        writeCost(objectives)
        writeLogHeader(objectives, var_set)
        logOutput(objectives, var_set)
        fun_lib.logObjectives(var_set['__LOG_OBJECTIVES_PATH'], objectives, var_set['__SORT_COSTS_BY'])

    except FileNotFoundError as f:
        # the simulation did not evens started
        logCrash(var_set['__VALUE_LOG_DIRNAME'], 'Simulation output file %s not found, quitting.' % f.filename)
        # objectives = [fun_lib.ObjectiveVar('Simulation failed', value = 1, costFunctionType=fun_lib.CostFunctionType.DistanceFromZero)]
        # writeCost(objectives)
        # logOutput(objectives)
    except fun_lib.SimulationResultIncompleteError as e:
        # most probably the simulation crashed mid-simulation. Lets provide a penalty for that
        objectives = [fun_lib.ObjectiveVar('Simulation failed', value = 10, costFunctionType=fun_lib.CostFunctionType.DistanceFromZero)]

        # comment out if you want to halt the optim
        writeCost(objectives)

        logCrash(var_set['__VALUE_LOG_DIRNAME'], e)
        logOutput(objectives, var_set)

processDyMatFile()
print('All dan!')