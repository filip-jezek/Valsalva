# This cost function just calls all the other cfs
import fun_lib
import importlib
from typing import Iterable

DEFAULT_TARGETVARS_TAG = 'All_supine'


def importCostFunction(location):
    """Returns CF module, location is expected at sibling folder"""

    # working dir is in  'identification/combined/debug/' or sibling folder
    path = '..\\..\\' + location + '\\cost_function.py'
    spec = importlib.util.spec_from_file_location(
        'cost_function', path)
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)
    return cf

def countTotalCost(objectives : Iterable[fun_lib.ObjectiveVar]):
    active_obj = sum(1 for o in objectives if o.costFunctionType is not fun_lib.CostFunctionType.Ignore)
    total_cost = sum(o.cost() for o in objectives)
    # weighing cost by numbr of objectives
    return total_cost / active_obj
    
def mapVarSet(vars_set, mapping):
    """ Creates new varset with mapped input vars_set from simulation to that used in cost function (aka original model)
    All other (non-mapped) stay the same
    """

    for key in vars_set.keys():
        if key in mapping:
            new_vars_set[key] = vars_set[mapping[key]]
        else:
            new_vars_set[key] = vars_set[key]

    return new_vars_set

def getObjectives(vars_set, targetsFileName = r'../targetValues_' + DEFAULT_TARGETVARS_TAG + '.txt'):
    """ Gets all objectives for our combined CF.
    Trials:
    
    valsalva sit
    # valsalva supine

    tilt 60
        HR
        brachial pressures
        LV volumes (SV and EF) 
    max exercise"""

    objectives = list()

    def buildCostObjective(name, cost_func_folder, mapping):
        cf = importCostFunction('valsalva')
        mapped_vars = mapVarSet(vars_set, mapping)
        objectives = cf.getObjectives(mapped_vars)
        cost = countTotalCost(objectives)
        costObjective = fun_lib.ObjectiveVar(name, value= cost, costFunctionType=fun_lib.CostFunctionType.DistanceFromZero)
        return costObjective

    # valsalva supine
    objectives.append(buildCostObjective('Valsalva sup', 'valsalva', {'model_var': 'costFunc_var'})
    objectives.append(buildCostObjective('Baseline', 'optimizeBaseline', {'model_var': 'costFunc_var'})
    objectives.append(buildCostObjective('Exercise', 'MaxExercise', {'model_var': 'costFunc_var'})
    
    return objectives