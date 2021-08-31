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
    vars_set['__draw_plots'] = False
    
    fun_lib.checkSimulationLength(vars_set['time'][-1],20, penalty=1000)
    # we have 3 intervals
    # baseline resting volumes and CO
    # baseline resting pressures and CO
    # exercise pressures and CO

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

    steady1 = 30
    steady2 = 45
    steady3 = 60
    interval = fun_lib.findInterval(time[0] + steady1-3, time[0] + steady1, time)

    # interval for averaging
    interval_ex = fun_lib.findInterval(time[0] + steady2-3, time[0] + steady2, time)
    # max exercise
    interval_mex = fun_lib.findInterval(time[0] + steady3-3, time[0] + steady3, time)
    # to observe how much the baseline fluctuates
    # steady_interval = fun_lib.findInterval(time[-1] - 30, time[-1] -5, time)

    # mitral valve flow ratio spontaneous:atrial contraction is about 2:1 
    # vla_1st_Peak_i = fun_lib.findInterval(39.8, 40, time)
    # vla_2nd_Peak_i = fun_lib.findInterval(39.2, 39.4, time)
    # vla_peak_frac = numpy.max(vars_set['q_mv'][vla_1st_Peak_i])/numpy.max(vars_set['q_mv'][vla_2nd_Peak_i])

    # interval for Q MV - must be a bit shorter to find further saddle after the peak
    # intervalQMV = fun_lib.findInterval(time[-1] - 2, time[-1] -1, time)
    # (q_mv_Ppassive, q_mv_Patrial) = fun_lib.calculateQ_MV(vars_set['q_mv'], time, intervalQMV)
    # (q_mv_Ppassive, q_mv_Patrial) = fun_lib.calculateQ_MV(vars_set['q_mv'], time, intervalQMV)
    # (E_A, A_s) = fun_lib.calculateEA(vars_set['q_mv'],vars_set['cardiac_cycle'], time, intervalQMV)
    # vla_peak_frac = q_mv_Ppassive/q_mv_Patrial
    # (atrial_kick, q_mv_saddle) = fun_lib.calculateQdot_mv(vars_set['V_LV'], vars_set['q_mv'], time, intervalQMV)
    # peak pressure drop on aortic valve
    # dp_av = numpy.max(vars_set['ascending_aorta'][interval]) - numpy.max(vars_set['P_LV'][interval])

    # Van Bortel 2012 siggest using 80 % of carotid to femoral distance
    # distance = vars_set['speedSegmentLength'][1]*0.8
    distance = 0.677*0.8
    # pwv = fun_lib.calculatePWV(
    #     vars_set['time'],
    #     vars_set['carotid_pressure'],
    #     vars_set['femoral_pressure'],
    #     distance)
    pwv = 10

    # build costs
    ov = [  
        # MRI
            ('ESV', numpy.min(vars_set['V_LV'][interval])/ml2SI, 51, None, 3),
            ('EDV', numpy.max(vars_set['V_LV'][interval])/ml2SI, 133, None, 3),
            ('CO', numpy.mean(vars_set['CO'][interval]) /lpm2SI, 5.87, None, 0.05),
            # RHC
            ('BPs', max(vars_set['brachial_pressure'][interval])/mmHg2SI, 153, None, 0.5),
            ('BPd', min(vars_set['brachial_pressure'][interval])/mmHg2SI, 83, None, 0.5),
            # ('BPk', numpy.mean(vars_set['renal_capillary'][interval]) /mmHg2SI, 20, None, 1, 1/mmHg2SI),
            ('Ppas', numpy.max(vars_set['P_pa'][interval])/mmHg2SI, 29, None, 5),
            ('Ppad', numpy.min(vars_set['P_pa'][interval])/mmHg2SI, 14, None, 5),
            ('Ppv', numpy.mean(vars_set['P_pv'][interval])/mmHg2SI, 13, None, 0.5),
            ('P_SRV', numpy.max(vars_set['heartComponent.ventricles.P_RV'][interval])/mmHg2SI, 32, None, 5),
            ('P_sv', numpy.min(vars_set['P_sv'][interval])/mmHg2SI, 4, None, 1),
            
            # 20W exercise
            ('BPs_Ex', max(vars_set['brachial_pressure'][interval_ex])/mmHg2SI, 159, None, 1),
            ('BPd_Ex', min(vars_set['brachial_pressure'][interval_ex])/mmHg2SI, 82, None, 2),
            ('Ppas_Ex', numpy.max(vars_set['P_pa'][interval_ex])/mmHg2SI, 42, None, 5),
            ('Ppad_Ex', numpy.min(vars_set['P_pa'][interval_ex])/mmHg2SI, 26, None, 5),
            ('Ppv_Ex', numpy.mean(vars_set['P_pv'][interval_ex])/mmHg2SI, 23, None, 1),
            ('SP_RV_Ex', numpy.max(vars_set['heartComponent.ventricles.P_RV'][interval_ex])/mmHg2SI, 45, None, 5),
            ('DP_RV_Ex', numpy.min(vars_set['P_sv'][interval_ex])/mmHg2SI, 6, None, 1),
            ('CO_Ex', numpy.mean(vars_set['CO'][interval_ex]) /lpm2SI, 11.5, None, 0.1),
          
            # Max exercise
            ('BPs_mEx', max(vars_set['brachial_pressure'][interval_mex])/mmHg2SI, 163, None, 2),
            ('BPd_mEx', min(vars_set['brachial_pressure'][interval_mex])/mmHg2SI, 113, None, 4),
            ('Ppas_mEx', numpy.max(vars_set['P_pa'][interval_mex])/mmHg2SI, 57, None, 5),
            ('Ppad_mEx', numpy.min(vars_set['P_pa'][interval_mex])/mmHg2SI, 37, None, 5),
            ('Ppv_mEx', numpy.mean(vars_set['P_pv'][interval_mex])/mmHg2SI, 25, None, 2),
            ('SP_RV_mEx', numpy.max(vars_set['heartComponent.ventricles.P_RV'][interval_mex])/mmHg2SI, 45, None, 0),
            ('DP_RV_mEx', numpy.min(vars_set['P_sv'][interval_mex])/mmHg2SI, 6, None, 0),
            ('CO_mEx', numpy.mean(vars_set['CO'][interval_mex]) /lpm2SI, 14.5, None, 0),

            # # general assumptions
            ('PWV', pwv, None, [8, 50], 0)            
            ]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], tolerance=o[4], k_p=1), ov))


    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives
