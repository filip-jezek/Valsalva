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
import sys
sys.path.append("..")
import fun_lib

# (prefix, cost_function folder, fmu name)
model_infos = [
        ('baseline', 'optimizeBaselineTriSegLumens', 'ADAN_0main_SystemicTree_Identification_Results_OlufsenTriSeg_0optimized_0steadyState_0init.fmu'), 
        ('exercise', 'MaxExercise', 'ADAN_0main_SystemicTree_Exercise_OlufsenTriseg_0Exercise.fmu'),
        ('tilt', 'optimizeTilt', 'ADAN_0main_0SystemicTree_0Tilt_0OlufsenTriSeg_0tiltable.fmu'),
        # ('tilt_d', 'optimizeTilt', 'ADAN_0main_0SystemicTree_0Tilt_0OlufsenTriSeg_0tiltable_dymsolv.fmu'),
        ('valsalva', 'valsalva', 'ADAN_0main_SystemicTree_Valsalva_OlufsenTriSeg_0valsalva.fmu')]

schedules = 'random_search_schedules.csv'
niter = 500
DRAW_PLOTS = True
# output = ['brachial_pressure']

def initSchedules(params:dict):
    with open(schedules, 'w') as file:
        file.write('Run, cost, %s\n' % ', '.join(params.keys()))

def writeSchedule(params:dict, cost, iter):
    if not os.path.exists(schedules):
        initSchedules(params)

    with open(schedules, 'a') as file:
        line = '%03d, %.6f, %s\n' % (iter, cost, ', '.join('%.3e' % p for p in params.values()))
        file.write(line)

def readInitParams(initparams_path = 'init_params_default_vals.csv', continue_params_path = 'random_search_schedules.csv'):

    if os.path.exists(continue_params_path):
        # lets continue where we ended
        print('Continue from last run in ' + continue_params_path)
        with open(continue_params_path, 'r') as file:
            # Run, cost, param1, param2...
            param_names = file.readline().rstrip('\n').split(',')[2:]
            # last line should have the lowest costs - strip it off the \n and floatize
            param_strings = (file.readlines()[-1].rstrip('\n').split(',')[2:])
            param_vals = (float(v) for v in param_strings)
        
        params_def = dict(zip(param_names, param_vals))
        return params_def

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

def logger(s):
    print('[FMI] ' + s)

def initFmus() -> list:
    fmu_init = []
    for model_info in model_infos:
        (prefix, cf_folder, fmu) = model_info
        # extract the FMU to a temporary directory
        unzipdir = fmpy.extract(fmu)
        # read the model description
        model_description = fmpy.read_model_description(unzipdir)
        # instantiate the FMU
        fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description, 'CoSimulation')
        # fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description, 'CoSimulation', debug_logging=True, fmi_call_logger=logger, visible=True)
        
        # TODO this might be obsolete atm
        cf = fun_lib.importCostFunction(dir = '..\\' + cf_folder + '\\')
        fmu_init.append((fmu_instance, model_description, cf, prefix))
           
    return fmu_init

def convertResult(result, prefix = '', exclude = []):
    var_set = {}
    for name in result.dtype.names:
        if name in exclude:
            var_set[name] = exclude[name]
        else:
            new_name = '.'.join([prefix, name])
            var_set[new_name] = result[name]

    return var_set



def runSimulations(fmu_init, params, run):

    # reset the FMU instance instead of creating a new one
    combined_cost_func = fun_lib.importCostFunction(dir = '..\\Combined' + '\\')
    var_set = {}

    for (fmu_instance, fmu_desc, _, prefix) in fmu_init:
        
        # TODO do the paralelization
        result = fmpy.simulate_fmu(fmu_instance.unzipDirectory,
                                stop_time=60,
                                start_values=params,
                                model_description=fmu_desc,
                                fmu_instance=fmu_instance
                                )

        var_set.update(convertResult(result, prefix=prefix, exclude=params))
        fmu_instance.reset()

    if DRAW_PLOTS:
        var_set['__draw_plots'] = True
        var_set['__plot_title'] = "Run %i" % (run)
        var_set['__saveFig_path'] = "%sFitFig_%03d.png" % (fun_lib.getSafeLogDir('//Schedules', safe_rollback=''), run)
        var_set['__showPlots'] = False

    
    objectives = combined_cost_func.getObjectives(var_set, targetsFolder = r"../../data/Valsalva/", top_level = '..\\..\\')

    return objectives

def runOptimization():

    fmu_init = initFmus()
    params = readInitParams()
    cost = 1e300
    # shuffledParams = copy.

    for iter in range(niter):
    # todo paralelize this

        cur_params = shuffleParams(params)
        
        try:
            objectives = runSimulations(fmu_init, cur_params, iter)
            cur_cost = fun_lib.countTotalSumCost(objectives)
            if cur_cost < cost:
                # woohoo, better costs!
                cost = cur_cost
                params = cur_params
                writeSchedule(params, cost, iter)
        except:
            writeSchedule(cur_params, 1.0, iter)

            
        
        


    # clean up
    for (fmu_instance, _, _, _) in fmu_init:
        # free the FMU instance and unload the shared library
        fmu_instance.freeInstance()
        # delete the temporary directory
        shutil.rmtree(fmu_instance.unzipDirectory, ignore_errors=True)

    print('that is all, folks')
    pass

if __name__ is '__main__':
    runOptimization()
