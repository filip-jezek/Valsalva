# import scipy.io as skipy

import fun_lib
# from statistics import mean
from matplotlib import pyplot as plt
import numpy

import math
# import DyMat
# from matplotlib import pyplot as plt

def buildObjectives(o, baseline, time):
    # (BP,      phase1, numpy.max   , 'ph1_peak'    ,  1/baseline, count),
    return None


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
    HR = vars_set['heart_rate']

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
    valsalva_end = fun_lib.getValsalvaEnd(valsalva_start, time, vars_set['thoracic_pressure'])
    
    # divide valsalva phases    
    # pre-valsalva
    phase0 = (0, valsalva_start)
    # overshoot phase at start of the valsalva
    phase1 = (valsalva_start, valsalva_start + 5)
    # min mean pressure during valsalva
    phase2 = (valsalva_start+2, valsalva_end)
    # drop and overshoot after valsalva release
    phase4 = (valsalva_end, valsalva_end + 7)
    # recovery - all is getting to normal
    phase5 = (time[-1] - 5, time[-1])


    # baseline is start to just before Valsalva
    baseline_interval = fun_lib.findInterval(phase0[0], phase0[1], time)
    baseline_bp = bp_mean[baseline_interval].mean()
    baseline_hr = HR[baseline_interval].mean()
    
    IGNORE = fun_lib.CostFunctionType.Ignore
    COUNT = fun_lib.CostFunctionType.Quadratic
    # construct objectives
    objectives = list()
    objectives.append(fun_lib.ObjectiveVar('baseline_bp', baseline_bp, costFunctionType=IGNORE))
    objectives.append(fun_lib.ObjectiveVar('baseline_hr', baseline_hr, costFunctionType=IGNORE))

    phase_values = [(BP     , phase1, (-1, 0), numpy.max  , baseline_bp, 'ph1_peak'    , COUNT),
                    (bp_mean, phase2, (0, -2), numpy.min  , baseline_bp, 'ph2_mean_min', COUNT),
                    (bp_mean, phase4,(-5, -5), numpy.max  , baseline_bp, 'ph2_max'     , COUNT),
                    (bp_mean, phase4, (0, -2), numpy.min  , baseline_bp, 'ph4_drop'    , COUNT),
                    (bp_mean, phase4, (2, 5) , numpy.max  , baseline_bp, 'ph4_ovrshoot', COUNT),
                    (bp_mean, phase5, 0      , numpy.mean , baseline_bp, 'ph5_recovery', COUNT),
                    (HR     , phase1, 0      , numpy.min  , baseline_hr, 'ph1_hr_min' , COUNT),
                    (HR     , phase4, (0, 0) , numpy.max  , baseline_hr, 'ph4_hr_max' , COUNT),
                    (HR     , phase4, (0, 3) , numpy.min  , baseline_hr, 'ph4_hr_drop', COUNT),
                    (HR     , phase5, 0      , numpy.mean , baseline_hr, 'ph5_hr_recovery', COUNT)                  ]

    time_values = [ (BP,      phase1, (-1, 0), numpy.argmax, 't_ph1_peak'    , COUNT),
                    (bp_mean, phase2, (0, -2), numpy.argmin, 't_ph2_mean_min', COUNT),
                    (bp_mean, phase4,(-5, -5), numpy.argmax, 't_ph2_max'     , COUNT),
                    (bp_mean, phase4, (0, -2), numpy.argmin, 't_ph4_drop'    , COUNT),
                    (bp_mean, phase4, (2, 5) , numpy.argmax, 't_ph4_ovrshoot', COUNT),
                    (HR     , phase1, 0      , numpy.argmin, 't_ph1_hr_min'  , COUNT),
                    (HR     , phase4, (0, 0) , numpy.argmax, 't_ph4_hr_max'  , COUNT),
                    (HR     , phase4, (0, 3) , numpy.argmin, 't_ph4_hr_drop' , COUNT)    ]

    def getInterval(phase, phase_offset, time):
        if phase_offset is None:
            offset = (0, 0)
        elif isinstance(phase_offset, int):
            offset = (phase_offset, phase_offset)
        else:
            offset = (phase_offset[0], phase_offset[1])
        
        return fun_lib.findInterval(phase[0] + offset[0], phase[1] + offset[1], time), offset

    def buildValueObjective(o):
        (sig, phase, phase_offset, fun, baseline, name, include_in_cost) = o
        interval, _ = getInterval(phase, phase_offset, time)
        value = fun(sig[interval])/baseline
        return fun_lib.ObjectiveVar(name, value=value, costFunctionType=include_in_cost)

    def buildTimeObjective(to):
        (sig, phase, phase_offset, fun, name, include_in_cost) = to
        interval, offset = getInterval(phase, phase_offset, time)        
        value = time[fun(sig[interval])] + offset[0]
        return fun_lib.ObjectiveVar(name, value=value, costFunctionType=include_in_cost)        

    # map the inputs to ObjectiveVar
    objectives.extend(map(buildValueObjective, phase_values))
    objectives.extend(map(buildTimeObjective, time_values))

      
    if '__draw_plots' in vars_set:
        plt.plot(time, BP, 'b')
        plt.plot(time, bp_mean, 'm')
        plt.plot(time, vars_set['thoracic_pressure'], 'g')
        plt.plot(time, vars_set['heart_rate'])

        # get objective by name shortcuts
        getObj = lambda name: fun_lib.getObjectiveByName(objectives, name).value
        
        plt.plot(phase0, [baseline_bp]*2, 'k')
        plt.plot(phase5, [getObj('ph5_recovery')*baseline_bp]*2, 'k')

        plt.plot(phase1, [getObj('ph1_peak'    )*baseline_bp]*2, 'k')
        plt.plot(phase2, [getObj('ph2_mean_min')*baseline_bp]*2, 'k')
        plt.plot(phase4, [getObj('ph4_drop'    )*baseline_bp]*2, 'k')

        plt.plot(getObj('t_ph1_peak'    ) + phase1[0], getObj('ph1_peak'    )*baseline_bp, '*r')
        plt.plot(getObj('t_ph2_mean_min') + phase2[0], getObj('ph2_mean_min')*baseline_bp, '*r')
        plt.plot(getObj('t_ph2_max'     ) + phase4[0], getObj('ph2_max'     )*baseline_bp, '*r')        
        plt.plot(getObj('t_ph4_drop'    ) + phase4[0], getObj('ph4_drop'    )*baseline_bp, '*r')
        plt.plot(getObj('t_ph4_ovrshoot') + phase4[0], getObj('ph4_ovrshoot')*baseline_bp, '*r')

        plt.plot(phase0, [baseline_hr]*2, 'c')
        plt.plot(getObj('t_ph1_hr_min' ) + phase1[0], getObj('ph1_hr_min' )*baseline_hr, '*m')
        plt.plot(getObj('t_ph4_hr_max' ) + phase4[0], getObj('ph4_hr_max' )*baseline_hr, '*m')
        plt.plot(getObj('t_ph4_hr_drop') + phase4[0], getObj('ph4_hr_drop')*baseline_hr, '*m')
        plt.plot(phase5, [getObj('ph5_hr_recovery')*baseline_hr]*2, 'c')


    fun_lib.UpdateObjectivesByValuesFromFile('targetValues.txt', objectives)

    return objectives
