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
    path = '..\\..\\' + location + '\\cost_function.py'
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


def getObjectives(vars_set:dict, targetsFileName = r'../targetValues_' + DEFAULT_TARGETVARS_TAG + '.txt'):
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
    objectives = list()
    dpi = 100
    fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = False, figsize = [1920/dpi, 1080/dpi], dpi=dpi)
    
    
    if '__plot_title' in vars_set:
        fig_title = vars_set['__plot_title']
        # get rid of the key
        del vars_set['__plot_title']
    else:
        fig_title = 'Debug run'

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

    def buildCostObjective(name, cost_func_folder):
        """
        the name is prefix in the large model as well
        """
        nonlocal axes_num
        nonlocal ax

        vars_set['__plot_axes'] = ax[axes_num]
        axes_num = axes_num + 1
        cf = importCostFunction(cost_func_folder)
        # mapped_vars = mapVarSet(vars_set, mapping)
        # filtering only vars relevant for that cost function, e.g. valsalva.brachial_pressure gives brachial_pressure
        
        mapped_vars = filterVarSet(vars_set, name + '.')
        objectives = cf.getObjectives(mapped_vars)

        cost = fun_lib.countTotalWeightedCost(objectives)
        costObjective = ObjectiveVar(name, value=cost, costFunctionType=CostFunctionType.DistanceFromZero)

        return costObjective

    objectives.append(buildCostObjective('baseline', 'optimizeBaseline'))
    objectives.append(buildCostObjective('tilt', 'optimizeTilt'))
    objectives.append(buildCostObjective('exercise', 'MaxExercise'))
    # open the data folder
    vars_set['__targetValuesFilename'] = r"../../../data/Valsalva/targetValues_All_supine.txt"
    objectives.append(buildCostObjective('valsalva', 'valsalva'))
    
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
    
    fig.suptitle(fig_title)

    pic_path = fun_lib.getSafeLogDir(r'..\Schedules') + fig_title
    plt.savefig(pic_path)
    plt.show()


                
    
    return objectives