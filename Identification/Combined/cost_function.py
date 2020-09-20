# This cost function just calls all the other cfs
import fun_lib
from fun_lib import ObjectiveVar, CostFunctionType
import importlib
from typing import Iterable
import matplotlib.pyplot as plt

DEFAULT_TARGETVARS_TAG = 'All_supine'


def importCostFunction(location):
    """Returns CF module, location is expected at sibling folder"""

    # working dir is in  'identification/combined/debug/' or sibling folder
    path = location + '\\cost_function.py'
    spec = importlib.util.spec_from_file_location(
        'cost_function', path)
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)
    return cf

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


def getObjectives(vars_set:dict, targetsFolder = r"../../../data/Valsalva/", top_level = '..\\..\\..\\'):
    """ Gets all objectives for our combined CF.
    Trials:
    
    valsalva sit
    # valsalva supine

    tilt 60
        HR
        brachial pressures
        LV volumes (SV and EF) 
    max exercise"""

    plt.close('all')
    dpi = 100
    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = False, figsize = [1920/dpi, 1080/dpi], dpi=dpi)
    
    
    if '__plot_title' in vars_set:
        fig_title = vars_set['__plot_title']
        # get rid of the key
        del vars_set['__plot_title']
    else:
        fig_title = 'Debug run'

    if '__saveFig_path' in vars_set and vars_set['__saveFig_path'] is not None:
        saveFig_path = vars_set['__saveFig_path']
        # prevent child cost functions to save figures
        del vars_set['__saveFig_path']
    else:
        saveFig_path = None
    
    if '__objectivesLog_path' in vars_set:
        objectivesLog_path = vars_set['__objectivesLog_path']
    else:
        objectivesLog_path = None

    def flat2gen(alist):
        # https://stackoverflow.com/questions/3172930/flattening-mixed-lists-in-python-containing-iterables-and-noniterables
        for item in alist:
            if isinstance(item, Iterable):
                for subitem in item: 
                    yield subitem
            else:
                yield item
            
    ax = list(flat2gen(axs))
    axes_num = 0

    cost_objectives = []
    all_objectives = []

    def buildCostObjective(name, cost_func_folder, weight):
        """
        the name is prefix in the large model as well
        """
        nonlocal axes_num
        nonlocal ax

        vars_set['__plot_axes'] = ax[axes_num]
        axes_num = axes_num + 1
        cf = importCostFunction(top_level + 'Identification\\' + cost_func_folder)
        # mapped_vars = mapVarSet(vars_set, mapping)
        # filtering only vars relevant for that cost function, e.g. valsalva.brachial_pressure gives brachial_pressure
        
        mapped_vars = filterVarSet(vars_set, name + '.')
        objectives = cf.getObjectives(mapped_vars)
        
        # normalize to variance type cost function
        for o in objectives:
            fun_lib.unifyCostFunc(o)

        cost = fun_lib.countTotalWeightedCost(objectives)
        costObjective = ObjectiveVar(name, value=cost, costFunctionType=CostFunctionType.DistanceFromZero, weight=weight)
        # add it to the set
        cost_objectives.append(costObjective)

        
        # prepare to compare to all objectives in a log
        def sumObjectives(o:ObjectiveVar):
            o.name = name + '.' + o.name
            o.weight = weight * o.weight
            return o

        all_objectives.extend(map(sumObjectives, objectives))

    # they append inside
    buildCostObjective('baseline', 'optimizeBaselineTriSegLumens', 10)
    buildCostObjective('tilt', 'optimizeTilt', 1)
    buildCostObjective('exercise', 'MaxExercise', 1)
    # open the data folder
    vars_set['__targetValuesFilename'] = targetsFolder + 'targetValues_All_supine.txt'
    buildCostObjective('valsalva', 'valsalva', 1)
    
    # def plotObjectives():
    #     fignums = plt.get_fignums()
    #     axes = []

    #     for f in fignums:
    #         if len(plt.figure(f).axes) > 0:
    #             axes.append(plt.figure(f).axes)
            
        
    #     # # unpack axes
    #     # axes = list(f.axes for f in figs if f is not None and f.axes is not None and len(f.axes) > 0)
    #     plt.close('all')
    #     plt.figure()

    #     _, ax = plt.subplots(nrows = len(axes), sharex = True)

    #     ax.append(axes)
    #     ax.

    # plotObjectives()
    if objectivesLog_path is not None:
        with open(objectivesLog_path, 'w') as file:
            total_cost = fun_lib.countTotalSumCost(all_objectives)
            file.write('Name, value, target, %%, total = %.6e\n' % total_cost)
            for o in sorted(all_objectives, key = lambda o: o.cost(), reverse=True):
                s = '%s,%.3e,%s,%02d\n' % (o.name, o.value, o.target_log(), round(o.cost()/total_cost*100))
                file.write(s)



    
    fig.suptitle('%s costs %.6f' % (fig_title, fun_lib.countTotalSumCost(cost_objectives)))

    if saveFig_path is not None:
        plt.savefig(saveFig_path)
    
    if '__showPlots' in vars_set and vars_set['__showPlots'] :
    # and fun_lib.getRunNumber() == 0:
        # zero means debug
        plt.show()


                
    
    return cost_objectives