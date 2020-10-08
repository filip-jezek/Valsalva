# test FMPy optim
import fmpy
import os
import random
from typing import Iterable
from collections import OrderedDict
import copy
import math
import shutil

# the fun_lib is up a level
sys.path.append("..")
import fun_lib


fmu = '..\\ADAN_0main_SystemicTree_Baseline_OlufsenTriSeg_0base.fmu'
schedules = 'random_search_schedules.csv'
niter = 500
DRAW_PLOTS = True
# output = ['brachial_pressure']

def initSchedules(params:dict):
    with open(schedules, 'w+') as file:
        file.write('Run, cost, %s\n' % ', '.join(params.keys()))

def writeSchedule(params:dict, cost, iter):
    if not os.path.exists(schedules):
        initSchedules(params)

    with open(schedules, 'w+') as file:
        line = '%03d, %.6f, %s\n' % (iter, cost, ', '.join('%.3e' % p for p in params.values()))
        file.write(line)

def readInitParams(initparams_path = 'init_params_default_vals.csv'):

    if not os.path.exists(initparams_path):
        raise FileNotFoundError('init_params_default_vals.csv not found. Generate that using CodeGen/manipulate_dsin.py')

    params_def = OrderedDict()
    with open(initparams_path) as file:
        for line in file.readlines():
            cols = line.split(',')
            params_def[cols[0]] = float(cols[1].rstrip('\n'))

    return params_def

def shuffleParams(params_def : dict, step = 0.01):

    shuffled_params = {}
    random.seed = 1
    for param_key, param_val in params_def.items():
        # shuffled_params[param_key] = random.uniform(param_val*(1-step), param_val*(1+step))
        shuffled_params[param_key] = random.normalvariate(mu = param_val, sigma = param_val*step)
    
    return shuffled_params


def runSimulation():
    # extract the FMU to a temporary directory
    unzipdir = fmpy.extract(fmu)
    # read the model description
    model_description = fmpy.read_model_description(unzipdir)
    # instantiate the FMU
    fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description, 'CoSimulation')

    cf = fun_lib.importCostFunction(dir = 'Combined\\')

    params = readInitParams()
    cost = 1e300
    # shuffledParams = copy.

    for iter in range(niter):
    # todo paralelize this

        
        cur_params = shuffleParams(params)

        # reset the FMU instance instead of creating a new one
        fmu_instance.reset()
        
        result = fmpy.simulate_fmu(unzipdir,
                                    stop_time=1,
                                    start_values=cur_params,
                                    model_description=model_description,
                                    fmu_instance=fmu_instance,
                                    fmi_call_logger=lambda s: print('[FMI] ' + s) )

        var_set = {}
        for name in result.dtype.names:
            var_set[name] = result[name]
            
        if DRAW_PLOTS:
            var_set['__draw_plots'] = True
            var_set['__plot_title'] = "Run %i" % (fun_lib.getRunNumber())
            var_set['__saveFig_path'] = "%sFitFig_%03d.png" % (fun_lib.getSafeLogDir('Schedules'), fun_lib.getRunNumber())

        objectives = cf.getObjectives(var_set, targetsFolder = r"../data/Valsalva/")
        cur_cost = fun_lib.countTotalSumCost(objectives)

        if cur_cost < cost:
            # woohoo, better costs!
            cost = cur_cost
            params = cur_params
            writeSchedule(params, cost, iter)

        print(result)
    

    # free the FMU instance and unload the shared library
    fmu_instance.freeInstance()

    # delete the temporary directory
    shutil.rmtree(unzipdir, ignore_errors=True)

    print('that is all, folks')
    pass

runSimulation()