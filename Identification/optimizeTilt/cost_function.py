import scipy.io as skipy
import fun_lib
from statistics import mean
# import DyMat
from matplotlib import pyplot as plt

# // calculate the cost function


mmHg2SI = 133.32
ml2SI = 1e-6
lpm2SI = 1e-3/60
bpm2SI = 1/60


def plotObjectives(vars_set, interval, objectives):
    ax = fun_lib.getAxes(vars_set)

    ax.plot(vars_set['time'], vars_set['brachial_pressure']/133.32, label='Brachial pressure mmHg')
    # ax.plot(vars_set['time'], vars_set['SV']*1e6, label='Stroke volume ml')    
    ax.plot(vars_set['time'], vars_set['CO']/lpm2SI, label='CO')   
    ax.plot(vars_set['time'], vars_set['HR']/bpm2SI, label='HR')    

    pack = (objectives, vars_set['time'], ax, interval)
    fun_lib.plotObjectiveTarget(pack, 'BPs', 1/133.32)
    fun_lib.plotObjectiveTarget(pack, 'BPd', 1/133.32)
    fun_lib.plotObjectiveTarget(pack, 'CO', 1/lpm2SI)
    fun_lib.plotObjectiveTarget(pack, 'HR', 1/bpm2SI)

    total_costs = fun_lib.countTotalWeightedCost(objectives)
    ax.set_title('Tilt costs %.6f' % total_costs)
    ax.set_ylim([0, 140])


def getObjectives(vars_set):

    fun_lib.checkSimulationLength(vars_set['time'][-1],110)

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89

    BPs_target = 113*mmHg2SI
    BPd_target = 90*mmHg2SI
    SV_target = 63*ml2SI
    HR_target = 81*bpm2SI
    CO_target = 5.1*lpm2SI
    phi_target = 0.3739
    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(t-5, t, vars_set['time'])

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval])/mmHg2SI, 122.4, None, 5),
            ('BPd', min(vars_set['brachial_pressure'][interval])/mmHg2SI, 86.4, None, 5),
            ('CO', fun_lib.avg(vars_set['CO'], interval)/lpm2SI, 4.608, None, .1),
            ('HR', fun_lib.avg(vars_set['HR'], interval)/bpm2SI, 78.12, None, 2)]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], tolerance = o[4]), ov))

    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives
