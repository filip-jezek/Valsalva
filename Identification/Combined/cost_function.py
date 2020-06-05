# This cost function just calls all the other cfs
import fun_lib
from fun_lib import ObjectiveVar, CostFunctionType
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

def countTotalCost(objectives : Iterable[ObjectiveVar]):
    active_obj = sum(1 for o in objectives if o.costFunctionType is not CostFunctionType.Ignore)
    
    def unifyCostFunc(o:ObjectiveVar):
        if o.costFunctionType == CostFunctionType.Linear:
            o.costFunctionType = CostFunctionType.LinearVariance
        elif o.costFunctionType == CostFunctionType.Quadratic:
            o.costFunctionType = CostFunctionType.QuadraticVariance
        # ignored and distanceFromZero are ignored
    
    # to have comparable cost function values one must have the stds ready
    map(unifyCostFunc, objectives)

    total_cost = sum(o.cost() for o in objectives)
    # weighing cost by numbr of objectives
    return total_cost / active_obj
    
def mapVarSet(vars_set, mapping):
    """ Creates new varset with mapped input vars_set from simulation to that used in cost function (aka original model)
    All other (non-mapped) stay the same
    """
    new_vars_set = {}
    for key in vars_set.keys():
        if key in mapping:
            new_vars_set[key] = vars_set[mapping[key]]
        else:
            new_vars_set[key] = vars_set[key]

    return new_vars_set

def filterVarSet(vars_set, filter_string):
    
    filtered_vars = {k[len(filter_string):]:v for (k,v) in vars_set.items() if k.startswith(filter_string)}
    # keep the control vars, just in case
    control_vars = {k:v for (k,v) in vars_set.items() if k.startswith("__") or k is "time"}
    filtered_vars.update(control_vars)
    
    return filtered_vars

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

    def buildCostObjective(name, cost_func_folder):
        """
        the name is prefix in the large model as well
        """
        cf = importCostFunction(cost_func_folder)
        # mapped_vars = mapVarSet(vars_set, mapping)
        # filtering only vars relevant for that cost function, e.g. valsalva.brachial_pressure gives brachial_pressure
        
        mapped_vars = filterVarSet(vars_set, name + '.')
        objectives = cf.getObjectives(mapped_vars)
        cost = countTotalCost(objectives)
        costObjective = ObjectiveVar(name, value= cost, costFunctionType=CostFunctionType.DistanceFromZero)
        return costObjective


    objectives.append(buildCostObjective('baseline', 'optimizeBaseline'))
    objectives.append(buildCostObjective('exercise', 'MaxExercise'))
    # open the data folder
    vars_set['__targetValuesFilename'] = r"../../../data/Valsalva/targetValues_All_sitting.txt"
    objectives.append(buildCostObjective('valsalva', 'valsalva'))
    objectives.append(buildCostObjective('tilt', 'optimizeTilt'))
    
    return objectives