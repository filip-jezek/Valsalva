import scipy.io as skipy
import fun_lib
from statistics import mean
import  matplotlib.pyplot as plt
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function

def plotObjectives(vars_set, interval, objectives):
    ax = fun_lib.getAxes(vars_set)

    ax.plot(vars_set['time'], vars_set['brachial_pressure']/133.32, label='Brachial pressure mmHg')
    ax.plot(vars_set['time'], vars_set['V_LV']*1000*1000, label='V_LV ml')
    ax.plot(vars_set['time'], vars_set['TEjection'], label='TEjection')
    ax.plot(vars_set['time'], vars_set['TFilling'], label='TFilling')

    pack = (objectives, vars_set, ax, interval)
    fun_lib.plotObjectiveTarget(pack, 'BPs', 1/133.32)
    fun_lib.plotObjectiveTarget(pack, 'EDV', 1e6)
    fun_lib.plotObjectiveTarget(pack, 'ESV', 1e6)
    fun_lib.plotObjectiveTarget(pack, 'Ts', 1, verticalalignment='top')
    fun_lib.plotObjectiveTarget(pack, 'Td', 1)
    
    total_costs = fun_lib.countTotalWeightedCost(objectives)
    ax.set_title('Exercise costs %.6f' % total_costs)


def getObjectives(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60

    BPs_target = 160*mmHg2SI
    EDV_target = 133*ml2SI
    ESV_target = 27*ml2SI

    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(t-5, t, vars_set['time'])

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval]), BPs_target, None, 1),
            ('EDV', max(vars_set['V_LV'][interval]), EDV_target, None, 1),
            ('ESV', min(vars_set['V_LV'][interval]), ESV_target, None, 1),
            ('Ts', max(vars_set['TEjection'][interval]), 0.166, [0.1, 0.2], 1e-4),
            ('Td', max(vars_set['TFilling'][interval]), 0.138, [0.1, 0.2], 1e-4),
        ]
    
    # make it a dict?
    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight = o[4]), ov))

    objectives[-1].costFunctionType = fun_lib.CostFunctionType.Linear
    objectives[-2].costFunctionType = fun_lib.CostFunctionType.Linear

    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
       plotObjectives(vars_set, interval, objectives)

    return objectives
