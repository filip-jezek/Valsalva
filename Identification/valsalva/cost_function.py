# import scipy.io as skipy

import fun_lib
# from statistics import mean
from matplotlib import pyplot as plt
import numpy

import math
# import DyMat
# from matplotlib import pyplot as plt

DEFAULT_TARGETVARS_TAG = 'All_supine'


def plotTargetValues(objectives, valsalva_start, valsalva_end, signal_end):
    # plot merged timecourse
    # valsalva_start = 20
    # valsalva_end = 35
    # signal_end = 55
    
    # All targets are relative to baseline, so we have to have that set first
    baseline_bp = fun_lib.getObjectiveByName(objectives, 'baseline_bp').value
    baseline_hr = fun_lib.getObjectiveByName(objectives, 'baseline_hr').value

    # just a shortcut
    def getTrgtVal(name):
        return fun_lib.getObjectiveByName(objectives, name).targetValue
    
    def getTrgtVar(name):
        return fun_lib.getObjectiveByName(objectives, name).std

    # color for BP and HR
    c_BP = 'g'
    c_HR = 'c'

    # baselines
    plt.errorbar((valsalva_start/2), getTrgtVal('baseline_bp'), yerr=getTrgtVar('baseline_bp'), fmt = c_BP, barsabove = True, zorder=3)
    plt.errorbar((valsalva_start/2), getTrgtVal('baseline_hr'), yerr=getTrgtVar('baseline_hr'), fmt = c_HR, barsabove = True, zorder=3)
    
    # recoveries
    plt.errorbar((signal_end - 2.5), getTrgtVal('ph5_recovery')*baseline_bp, yerr=getTrgtVar('ph5_recovery')*baseline_bp, fmt = c_BP, barsabove = True, zorder=3)
    plt.errorbar((signal_end - 2.5), getTrgtVal('ph5_hr_recovery')*baseline_hr, yerr=getTrgtVar('ph5_hr_recovery')*baseline_hr, fmt = c_HR, barsabove = True, zorder=3)

    def plotMetric(t_val, t_offset, val, baseline, color):
        val_mean = getTrgtVal(val)*baseline
        val_std = getTrgtVar(val)*baseline
        t_mean = getTrgtVal(t_val) + t_offset
        t_std = getTrgtVar(t_val) 
        # zorder 3 is a workaround to show errorbars above the plots
        plt.errorbar(t_mean, val_mean, yerr= val_std, xerr=t_std, fmt = color, barsabove = True, zorder=3, linewidth = 2)

    def plotBPMetric(t_val, t_offset, val, color = c_BP):
        plotMetric(t_val, t_offset, val, baseline_bp, color)

    def plotHRMetric(t_val, t_offset, val, color = c_HR):
        plotMetric(t_val, t_offset, val, baseline_hr, color)    

    plotBPMetric('t_ph1_peak', valsalva_start, 'ph1_peak')
    plotBPMetric('t_ph2_mean_min', valsalva_start, 'ph2_mean_min')
    plotBPMetric('t_ph2_max', valsalva_end, 'ph2_max')
    plotBPMetric('t_ph4_drop', valsalva_end, 'ph4_drop')
    plotBPMetric('t_ph4_ovrshoot', valsalva_end, 'ph4_ovrshoot')

    plotHRMetric('t_ph1_hr_min', valsalva_start, 'ph1_hr_min')
    plotHRMetric('t_ph4_hr_max', valsalva_end, 'ph4_hr_max')
    plotHRMetric('t_ph4_hr_drop', valsalva_end, 'ph4_hr_drop')
    pass


def getObjectives(vars_set, targetsFileName = r'../targetValues_' + DEFAULT_TARGETVARS_TAG + '.txt'):
    """ Returns dict of objectives
    control variables:
    __targetValuesFilename
    __plot_title
    __saveFig_path
    __draw_plots"""
    # some control variables are not present in case of identification
    if '__targetValuesFilename' not in vars_set or vars_set['__targetValuesFilename'] is None:
        vars_set['__targetValuesFilename'] = targetsFileName

    # if '__draw_plots' not in vars_set:
    #     vars_set['__draw_plots'] = True
    

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    # ml2SI = 1e-6
    # lpm2SI = 1e-3/60
    BPM2SI = 1/60

    BP = vars_set['brachial_pressure']
    # HR = vars_set['heartRate.HR']
    TP = vars_set['thoracic_pressure']
    time = vars_set['time']

    # make sure its in non-SI units
    if numpy.mean(BP) > 200:
        # convert to SI units
        BP = BP/mmHg2SI
        # HR = HR/BPM2SI
        TP = TP/mmHg2SI

    dt = vars_set['time'][2] - vars_set['time'][1]
    assert dt > 0, "The simulation must not store values at events."

    
    peaks  = fun_lib.getPeaks(BP, dt)
    HR = fun_lib.getHR(BP, peaks, dt)/BPM2SI
    bp_mean = fun_lib.getMeanRR(BP, peaks)

    # find valsalva start and end
    valsalva_start = fun_lib.getValsalvaStart(time, TP, threshold=10)
    valsalva_end = fun_lib.getValsalvaEnd(valsalva_start, time, TP, threshold=10)
    valsalva = (valsalva_start, valsalva_end)
    
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
                    (bp_mean, phase4,(-7, -7), numpy.max  , baseline_bp, 'ph2_max'     , COUNT),
                    (bp_mean, phase4,(-2, -2), numpy.min  , baseline_bp, 'ph4_drop'    , COUNT),
                    (bp_mean, phase4, (2, 5) , numpy.max  , baseline_bp, 'ph4_ovrshoot', COUNT),
                    (bp_mean, phase5, 0      , numpy.mean , baseline_bp, 'ph5_recovery', COUNT),
                    (HR     , phase1, 0      , numpy.min  , baseline_hr, 'ph1_hr_min' , COUNT),
                    (HR     , phase4, (0, 0) , numpy.max  , baseline_hr, 'ph4_hr_max' , COUNT),
                    (HR     , phase4, (0, 3) , numpy.min  , baseline_hr, 'ph4_hr_drop', COUNT),
                    (HR     , phase5, 0      , numpy.mean , baseline_hr, 'ph5_hr_recovery', COUNT)                  ]

    time_values = [ (BP,      phase1, (-1, 0), numpy.argmax, 't_ph1_peak'    , COUNT),
                    (bp_mean, phase2, (0, -2), numpy.argmin, 't_ph2_mean_min', COUNT),
                    (bp_mean, phase4,(-7, -7), numpy.argmax, 't_ph2_max'     , COUNT),
                    (bp_mean, phase4, (-2, -2), numpy.argmin, 't_ph4_drop'    , COUNT),
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

    if '__targetValuesFilename' in vars_set and vars_set['__targetValuesFilename'] is not None:
        fun_lib.updateObjectivesByValuesFromFile(vars_set['__targetValuesFilename'], objectives)
      
    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plt.figure()

        if '__plot_title' in vars_set:
            plt.title(vars_set['__plot_title'])

        plt.plot(time, BP, 'b')
        plt.plot(time, bp_mean, 'm')
        plt.plot(time, TP, 'g')
        plt.plot(time, HR)
        # plt.plot(time, vars_set['heartRate.HR'], '--')


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

        plotTargetValues(objectives, valsalva_start, valsalva_end, time[-1])
        # plt.show(block = False)

        if '__saveFig_path' in vars_set and vars_set['__saveFig_path'] is not None:
            plt.savefig(vars_set['__saveFig_path'], dpi = 300)
        

    return objectives
