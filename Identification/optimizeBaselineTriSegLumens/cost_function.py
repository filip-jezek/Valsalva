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

    
    fun_lib.checkSimulationLength(vars_set['time'][-1],10, penalty=1000)

    # def 'johan' 


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

    # interval for averaging
    interval = fun_lib.findInterval(time[-1] - 2, time[-1], time)
    # to observe how much the baseline fluctuates
    steady_interval = fun_lib.findInterval(time[-1] - 30, time[-1] -5, time)

    # mitral valve flow ratio spontaneous:atrial contraction is about 2:1 
    # vla_1st_Peak_i = fun_lib.findInterval(39.8, 40, time)
    # vla_2nd_Peak_i = fun_lib.findInterval(39.2, 39.4, time)
    # vla_peak_frac = numpy.max(vars_set['q_mv'][vla_1st_Peak_i])/numpy.max(vars_set['q_mv'][vla_2nd_Peak_i])

    # interval for Q MV - must be a bit shorter to find further saddle after the peak
    intervalQMV = fun_lib.findInterval(time[-1] - 2, time[-1] -1, time)
    # (q_mv_Ppassive, q_mv_Patrial) = fun_lib.calculateQ_MV(vars_set['q_mv'], time, intervalQMV)
    # (q_mv_Ppassive, q_mv_Patrial) = fun_lib.calculateQ_MV(vars_set['q_mv'], time, intervalQMV)
    (q_mv_Ppassive, q_mv_Patrial, ea) = fun_lib.calculateEA(vars_set['q_mv'],vars_set['cardiac_cycle'], time, intervalQMV)
    # vla_peak_frac = q_mv_Ppassive/q_mv_Patrial
    # (atrial_kick, q_mv_saddle) = fun_lib.calculateQdot_mv(vars_set['V_LV'], vars_set['q_mv'], time, intervalQMV)
    # peak pressure drop on aortic valve
    dp_av = numpy.max(vars_set['ascending_aorta'][interval]) - numpy.max(vars_set['P_LV'][interval])

    # Van Bortel 2012 siggest using 80 % of carotid to femoral distance
    # distance = vars_set['speedSegmentLength'][1]*0.8
    distance = 0.677*0.8
    pwv = fun_lib.calculatePWV(
        vars_set['time'],
        vars_set['carotid_pressure'],
        vars_set['femoral_pressure'],
        distance)

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval])/mmHg2SI, 120, None, 1),
            ('BPd', min(vars_set['brachial_pressure'][interval])/mmHg2SI, 80, None, 1),
            ('EDV', numpy.max(vars_set['V_LV'][interval])/ml2SI, 150, None, 3),
            ('ESV', numpy.min(vars_set['V_LV'][interval])/ml2SI, 60, None, 3),
            ('ESV_la', numpy.min(vars_set['V_la'][interval])/ml2SI, 12, None, 5),
            ('EDV_la', numpy.max(vars_set['V_la'][interval])/ml2SI, 41, None, 5),
            ('dp_av', dp_av/mmHg2SI, None, [-5, 5], 0.1),
            ('E/A', ea, 1.7, None, 0.5),
            # ('Q_MV_f', vla_peak_frac*1, None, [1.5, 2], 0),
            # ('Qdot_mv', atrial_kick*1, None, [0.2, 0.3], 0),
            # ('q_mv_sad', q_mv_saddle*100, None, [0, q_mv_Patrial*20], 0),
            # ('SL_max', max(vars_set['SLo_max'][interval]), 2.2, None, 0, 1),
            # ('SL_min', min(vars_set['SLo_min'][interval]), 1.75, None, 0, 1),            
            # ('HR', numpy.mean(vars_set['HR'][interval]) , 64*bpm2SI, None, 1, 1/bpm2SI), 
# set by assumption and loop closed
#            ('HR', numpy.mean(vars_set['HR'][interval]), HR_target, None, 1), 
# set by EDV and ESV
#            ('EF', fun_lib.calculateEF(vars_set['V_LV'][interval]), EF_target, None, 1),
# set by HR and EDV AND ESV, just emphasized here
            ('CO', sum(vars_set['CO'][interval]) / len(interval)/lpm2SI, 5.76, None, 0.1),
            # ('BPk', sum(vars_set['renal_capillary'][interval]) / len(interval)/mmHg2SI, 20, None, 1, 1/mmHg2SI),
            ('Ppas', numpy.max(vars_set['P_pa'][interval])/mmHg2SI, 20.5, None, 5),
            ('Ppad', numpy.min(vars_set['P_pa'][interval])/mmHg2SI, 8.8, None, 5),
            ('Ppv', numpy.mean(vars_set['P_pv'][interval])/mmHg2SI, 8, None, 2),
            # ('EDP', numpy.min(vars_set['P_LV'][interval]), None, [6*mmHg2SI, 12*mmHg2SI], 1e-3),
            # ('P_MV_o', numpy.mean(vars_set['P_MV_o'][interval])/mmHg2SI, 4.5*mmHg2SI, None, 0, 1/mmHg2SI),
            # ('P_MV_c', numpy.mean(vars_set['P_MV_c'][interval])/mmHg2SI, 6*mmHg2SI, None, 0, 1/mmHg2SI),
            # ('BPMeanStd', numpy.std(vars_set['brachial_pressure_mean'][steady_interval]), None, [0, 4*mmHg2SI], 0, 1/mmHg2SI),
            ('Ts', max(vars_set['TEjection'][interval]), 0.292, None, 0.04),
            ('PWV', pwv, None, [5, 10], 1)            ]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], tolerance=o[4], k_p=1), ov))


    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives
