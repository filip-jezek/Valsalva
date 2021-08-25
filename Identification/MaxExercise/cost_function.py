import scipy.io as skipy
import fun_lib
import numpy
from statistics import mean
import  matplotlib.pyplot as plt
# import DyMat
# from matplotlib import pyplot as plt

mmHg2SI = 133.32
ml2SI = 1e-6
lpm2SI = 1e-3/60

def plotObjectives(vars_set, interval, objectives):
    
    ax = fun_lib.getAxes(vars_set)

    ax.plot(vars_set['time'], vars_set['brachial_pressure']/133.32, label='Brachial pressure mmHg')
    # ax.plot(vars_set['time'], vars_set['V_LV']*1000*1000, label='V_LV ml')
    ax.plot(vars_set['time'], vars_set['CO']/lpm2SI, label='CO')

    # ax.plot(vars_set['time'], vars_set['TEjection'], label='TEjection')
    # ax.plot(vars_set['time'], vars_set['TFilling'], label='TFilling')

    pack = (objectives, vars_set['time'], ax, interval)
    fun_lib.plotObjectiveTarget(pack, 'BPs', 1/133.32)
    # fun_lib.plotObjectiveTarget(pack, 'EDV', 1e6)
    fun_lib.plotObjectiveTarget(pack, 'ESV', 1e6)
    fun_lib.plotObjectiveTarget(pack, 'Ts', 1, verticalalignment='top')
    fun_lib.plotObjectiveTarget(pack, 'Td', 1)
    fun_lib.plotObjectiveTarget(pack,'Ppa', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'Ppa_s', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'Ppa_d', 1/133.32)        
    fun_lib.plotObjectiveTarget(pack,'Ppv', 1/133.32)
    fun_lib.plotObjectiveLimit(pack, 'CO', 1/lpm2SI, 'lower')
    
    total_costs = fun_lib.countTotalWeightedCost(objectives)
    ax.set_title('Exercise costs %.6f' % total_costs)
    ax.set_ylim([0, 165])



def getObjectives(vars_set):

    fun_lib.checkSimulationLength(vars_set['time'][-1],30, penalty=1000)

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89

    # t = vars_set['time'][-1]
    t = 30
    interval = fun_lib.findInterval(t-3, t, vars_set['time'])

    BPs_target = 194.4*mmHg2SI
    BPd_target = 88*mmHg2SI
    EDV_target = 133*ml2SI
    ESV_target = 30*ml2SI
    CO_min = 15.5*lpm2SI
    # Ppa_target = 34*mmHg2SI # (Kovacs Eur Respir J 2009)
    Ppv_target = 12*mmHg2SI # (Kovacs Eur Respir J 2009)
    co = fun_lib.calculateAlternativeCO(vars_set, interval)
    Ppv_target =co*60*1e3*1.2 # Eisman 2018 (PMID 29695381) shows a relation CO to PCWP of 1.2 mmHg/(L/min) in healthy control
    Ppas_target = 34*mmHg2SI # (Kovacs Eur Respir J 2009)
    
    Ppad_target = 15*mmHg2SI # (Kovacs Eur Respir J 2009)

    tm = 60
    intervalm = fun_lib.findInterval(tm-3, tm, vars_set['time'])
    co_max = fun_lib.calculateAlternativeCO(vars_set, intervalm)
    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval])/mmHg2SI, 194.4, None, 3),
            ('BPd', min(vars_set['brachial_pressure'][interval])/mmHg2SI, 81.16, None, 10),
            # ('EDV', max(vars_set['V_LV'][interval]), EDV_target, None, 1),
            # ('ESV', min(vars_set['V_LV'][interval]), ESV_target, None, 1e-3),
            ('CO', co/lpm2SI, 16.17, None, 0.1),
            ('CO_max', co_max/lpm2SI, 20.6, None, 1),
            ('Ts', max(vars_set['TEjection'][interval]), 0.196, None, .05),
            # ('Td', max(vars_set['TFilling'][interval]), 0.18, [0.18*0.5, 0.18*1.5], 1e-4),
            ('EF', fun_lib.calculateEF(vars_set['V_LV'][interval]), 0.78, None, .1),
            # ('HR', numpy.mean(vars_set['HR'][interval])*60, 154, None, 0, 60),
            # ('Ppa', numpy.mean(vars_set['P_pa'][interval]), Ppa_targ)et, None, .1),
            # ('Ppv', numpy.mean(vars_set['P_pv'][interval]), Ppv_target, None, .1),
            ('Ppa_s', numpy.max(vars_set['P_pa'][interval]), Ppas_target, None, 10),
            # ('Ppa_d', numpy.min(vars_set['P_pa'][interval]), Ppad_target, None, 0.1, 1/mmHg2SI),
            ('Ppv', numpy.mean(vars_set['P_pv'][interval])/mmHg2SI, Ppv_target, None, 2),
            ('Pcap', numpy.mean(vars_set['exercised_capillary'][interval])/mmHg2SI, None, [10, 40], 10),
        ]
    
    # make it a dict?
    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], tolerance = o[4], k_p=1), ov))

    # objectives[-1].costFunctionType = fun_lib.CostFunctionType.Linear
    # objectives[-2].costFunctionType = fun_lib.CostFunctionType.Linear

    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
       plotObjectives(vars_set, interval, objectives)

    return objectives
