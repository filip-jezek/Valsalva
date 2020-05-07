# import scipy.io as skipy

import fun_lib
# from statistics import mean
from matplotlib import pyplot as plt

import math
# import DyMat
# from matplotlib import pyplot as plt

def getObjectives(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60

    BP = vars_set['brachial_pressure']
    dt = vars_set['time'][2] - vars_set['time'][1]
    assert dt > 0, "The simulation must not store values at events."

    bp_mean = fun_lib.getMeanRR(BP, dt)
    time = vars_set['time']

    # bp_mean1 = getMeanSavGol(BP, dt)
    # bp_mean2 = getMeanRR(BP, dt)
    # plt.plot(time, BP)
    # plt.plot(time, bp_mean1)
    # plt.plot(time[0:len(bp_mean2)], bp_mean2)
    # sigf = ss.medfilt(BP, int(1/dt)+1)
    # plt.plot(time[0:len(sigf)], sigf)

    # plt.legend(['bp', 'mean1', 'mean2', 'medfilt'])
    
    # plt.show()


    # find valsalva start and end
    valsalva_start = fun_lib.getValsalvaStart(time, vars_set['thoracic_pressure'])
    valsalva_end = fun_lib.getValsalvaEnd(valsalva_start + 5, time, vars_set['thoracic_pressure'])
    # baseline is start to just before Valsalva
    baseline = bp_mean[fun_lib.findInterval(0, valsalva_start, time)].mean()
    
    # systolic peak at start of valsalva
    ph1_peak = max(BP[fun_lib.findInterval(valsalva_start, valsalva_start + 5, time)])
    # min mean pressure during valsalva
    ph2_mean_min = min(bp_mean[fun_lib.findInterval(valsalva_start, valsalva_end-2, time)])
    # peak minimal after valsalva release
    ph3_mean_min = min(bp_mean[fun_lib.findInterval(valsalva_end, valsalva_end + 10, time)])
    # peak maximal after valsalva
    ph4_mean_max = max(BP[fun_lib.findInterval(valsalva_end, valsalva_end + 15, time)])
    # mean recovery
    ph5_recovery_mean = bp_mean[fun_lib.findInterval(valsalva_end+15, valsalva_end + 20, time)].mean()


    # construct objectives
    objectives = list()
    objectives.append(fun_lib.ObjectiveVar('baseline', baseline, costFunctionType=fun_lib.CostFunctionType.Ignore))

    # build costs, relative to baseline
    ov = [  ('ph1', ph1_peak/baseline, None),
            ('ph2', ph2_mean_min/baseline, None),
            ('ph3', ph3_mean_min/baseline, None),
            ('ph4', ph4_mean_max/baseline, None),
            ('ph5', ph5_recovery_mean/baseline, None),
            ]
    # map the inputs to ObjectiveVar
    objectives.extend(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2]), ov))
  
    if '__draw_plots' in vars_set:
        plt.plot(time, BP, 'b')
        plt.plot(time, bp_mean, 'r')
        plt.plot((valsalva_start - 10, valsalva_start), [baseline]*2, 'k')
        plt.plot((valsalva_start, valsalva_start + 5), [ph1_peak]*2, 'k')
        plt.plot((valsalva_start, valsalva_end - 2), [ph2_mean_min]*2, 'k')
        plt.plot((valsalva_end, valsalva_end+10), [ph3_mean_min]*2, 'k')
        plt.plot((valsalva_end, valsalva_end+15), [ph4_mean_max]*2, 'k')
        plt.plot((valsalva_end+15, valsalva_end+20), [ph5_recovery_mean]*2, 'k')
        plt.draw()

    return objectives
