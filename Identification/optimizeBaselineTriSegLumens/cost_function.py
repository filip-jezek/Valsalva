import scipy.io as skipy
import fun_lib
import re
import matplotlib.pyplot as plt
import numpy
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function
def plotObjectives(vars_set, interval, objectives):
    if '__plot_axes' in vars_set:
        ax = vars_set['__plot_axes']
    else:
        fig = plt.figure()
        ax = fig.subplots()

    ax.plot(vars_set['time'], vars_set['brachial_pressure']/133.32, label='Brachial pressure mmHg')
    ax.plot(vars_set['time'], vars_set['CO']*1000*60, label='CO l/min')
    ax.plot(vars_set['time'], vars_set['renal_capillary']/133.32, label='Capillary pressure')
    # ax.plot([vars_set['time'][interval[0]], vars_set['time'][interval[-1]]], [fun_lib.getObjectiveByName(objectives, 'EF').value*100]*2, label='EF')
    ax.plot([vars_set['time'][interval[0]], vars_set['time'][interval[-1]]], [fun_lib.getObjectiveByName(objectives, 'PWV').value*1]*2, label='PWV')

    # bounds
    pack = (objectives, vars_set['time'], ax, interval)
    fun_lib.plotObjectiveTarget(pack,'BPs', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'BPd', 1/133.32)
    # fun_lib.plotObjectiveTarget(pack,'CO', 1000*60)
    fun_lib.plotObjectiveTarget(pack,'BPk', 1/133.32)
    # fun_lib.plotObjectiveTarget(pack,'EF', 100)
    fun_lib.plotObjectiveTarget(pack,'HR', 60, verticalalignment='top') 
    # fun_lib.plotObjectiveTarget(pack,'Ppa', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'Ppas', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'Ppad', 1/133.32)        
    fun_lib.plotObjectiveTarget(pack,'Ppv', 1/133.32)
    fun_lib.plotObjectiveLimit(pack, 'PWV', 1, 'lower', verticalalignment='top')
    

    total_costs = fun_lib.countTotalWeightedCost(objectives)
    ax.set_title('Baseline costs %.6f' % total_costs)
    ax.set_ylim([0, 140])
    # ax.set_xlim(60, 120)



def getObjectives(vars_set):

    
    fun_lib.checkSimulationLength(vars_set['time'][-1],90)


    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60
    bpm2SI = 1/60

    # BPs_target = 120*mmHg2SI
    # BPd_target = 80*mmHg2SI
    # BPk_target = 20*mmHg2SI
    # CO_target = 5.76*lpm2SI
    # EF_target = 0.6
    # HR_target = 60*bpm2SI
    # Ppa_target = 14*mmHg2SI # (Kovacs Eur Respir J 2009)
    # Ppv_target = 8*mmHg2SI # (Kovacs Eur Respir J 2009)
    # Ppas_target = 20.8*mmHg2SI # (Kovacs Eur Respir J 2009)
    # Ppad_target = 8.8*mmHg2SI # (Kovacs Eur Respir J 2009)
    # pwv_bounds = [5, 10]

    time = vars_set['time']

    # interval = fun_lib.findInterval(time[0] + 43, time[0] + 48, time)
    interval = fun_lib.findInterval(time[-1] - 8, time[-1] - 3, time)
    # to observe how much the baseline fluctuates
    steady_interval = fun_lib.findInterval(time[-1] - 30, time[-1] -5, time)

    # mitral valve flow ratio spontaneous:atrial contraction is about 2:1 
    # vla_1st_Peak_i = fun_lib.findInterval(39.8, 40, time)
    # vla_2nd_Peak_i = fun_lib.findInterval(39.2, 39.4, time)
    # vla_peak_frac = numpy.max(vars_set['q_mv'][vla_1st_Peak_i])/numpy.max(vars_set['q_mv'][vla_2nd_Peak_i])
    (q_mv_Ppassive, q_mv_Patrial) = fun_lib.calculateQ_MV(vars_set['q_mv'], time, interval)
    vla_peak_frac = q_mv_Ppassive/q_mv_Patrial
    (atrial_kick, q_mv_saddle) = fun_lib.calculateQdot_mv(vars_set['V_LV'], vars_set['q_mv'], time, interval)

    # Van Bortel 2012 siggest using 80 % of carotid to femoral distance
    # distance = vars_set['speedSegmentLength'][1]*0.8
    distance = 0.677*0.8
    pwv = fun_lib.calculatePWV2(
        vars_set['time'],
        vars_set['carotid_pressure'],
        vars_set['femoral_pressure'],
        distance)

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval]), 120*mmHg2SI, None, 10),
            ('BPd', min(vars_set['brachial_pressure'][interval]), 80*mmHg2SI, None, 10),
            ('EDV', numpy.max(vars_set['V_LV'][interval]), 150*ml2SI, None, 1e-1),
            ('ESV', numpy.min(vars_set['V_LV'][interval]), 60*ml2SI, None, 1e-1),
            ('ESV_la', numpy.min(vars_set['V_la'][interval]), 12*ml2SI, None, 1e-2),
            ('EDV_la', numpy.max(vars_set['V_la'][interval]), 41*ml2SI, None, 1e-2),
            ('Q_MV_f', vla_peak_frac, None, [1.5, 2], 1e-3),
            ('Qdot_mv', atrial_kick, None, [0.2, 0.3], 1e-3),
            ('q_mv_sad', q_mv_saddle, None, [0, q_mv_Patrial*0.2], 1e-3),
            ('SL_max', max(vars_set['SLo_max'][interval]), 2.2, None, 0),
            ('SL_min', min(vars_set['SLo_min'][interval]), 1.75, None, 0),            
            ('HR', numpy.mean(vars_set['HR'][interval]) , 64*bpm2SI, None, 1), 
# set by assumption and loop closed
#            ('HR', numpy.mean(vars_set['HR'][interval]), HR_target, None, 1), 
# set by EDV and ESV
#            ('EF', fun_lib.calculateEF(vars_set['V_LV'][interval]), EF_target, None, 1),
# set by HR and EDV AND ESV, just emphasized here
            ('CO', sum(vars_set['CO'][interval]) / len(interval), 5.76*lpm2SI, None, 10),
            ('BPk', sum(vars_set['renal_capillary'][interval]) / len(interval), 20*mmHg2SI, None, .1),
            ('Ppas', numpy.max(vars_set['P_pa'][interval]), 20.5*mmHg2SI, None, .1),
            ('Ppad', numpy.min(vars_set['P_pa'][interval]), 8.8*mmHg2SI, None, .1),
            ('Ppv', numpy.mean(vars_set['P_pv'][interval]), 8*mmHg2SI, None, .1),
            ('EDP', numpy.min(vars_set['P_LV'][interval]), None, [6*mmHg2SI, 12*mmHg2SI], 1e-3),
            ('BPMeanStd', numpy.std(vars_set['brachial_pressure_mean'][steady_interval]), None, [0, 4*mmHg2SI], 1e-2),
            ('PWV', pwv, None, [5, 10], 1)            ]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight=o[4], k_p=1e3), ov))


    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives
