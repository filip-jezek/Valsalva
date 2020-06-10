import scipy.io as skipy
import fun_lib
from statistics import mean
# import DyMat
from matplotlib import pyplot as plt

# // calculate the cost function
def plotObjectives(vars_set, interval, objectives):
    ax = fun_lib.getAxes(vars_set)

    ax.plot(vars_set['time'], vars_set['brachial_pressure']/133.32, label='Brachial pressure mmHg')
    ax.plot(vars_set['time'], vars_set['SV']*1e6, label='Stroke volume ml')    
    ax.plot(vars_set['time'], vars_set['phi_baro']*100, label='phi %')    

    pack = (objectives, vars_set, ax, interval)
    fun_lib.plotObjectiveTarget(pack, 'BPs', 1/133.32)
    fun_lib.plotObjectiveTarget(pack, 'BPd', 1/133.32)
    fun_lib.plotObjectiveTarget(pack, 'SV', 1e6)
    fun_lib.plotObjectiveTarget(pack, 'phi', 100)

    total_costs = fun_lib.countTotalWeightedCost(objectives)
    ax.set_title('Tilt costs %.6f' % total_costs)

def getObjectives(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60

    BPs_target = 113*mmHg2SI
    BPd_target = 90*mmHg2SI
    SV_target = 63*ml2SI
    phi_target = 0.3739
    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(t-5, t, vars_set['time'])

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval]), BPs_target, None, 1),
            ('BPd', min(vars_set['brachial_pressure'][interval]), BPd_target, None, 0.1),
            ('SV', fun_lib.avg(vars_set['SV'], interval), SV_target, None, 0.1),
            ('phi', fun_lib.avg(vars_set['phi_baro'], interval), phi_target, None, 1)]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight = o[4]), ov))

    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives
