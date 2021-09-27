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
    
    fun_lib.checkSimulationLength(vars_set['time'][-1],15, penalty=1000)
    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60
    bpm2SI = 1/60

    time = vars_set['time']

    interval = fun_lib.findInterval(time[-1] -3, time[-1], time)
  
    
    # ratio of CO calculated from LV volumes and other methods.
    # Used to scale the LV volumes to reasonable values
    LV_upscale = 4470/(68-26)/64 

    # values from Jones et al. 2021, HFrEF patients with excluded outliers, elaving patients 4, 7 and 8
    # build costs
    ov = [  
            # ECHO
            ('EDV', numpy.max(vars_set['V_LV'][interval])/ml2SI, 129, None, 33),
            ('ESV', numpy.min(vars_set['V_LV'][interval])/ml2SI, 65, None, 27),
            ('CO', numpy.mean(vars_set['CO'][interval]) /lpm2SI, (4.9 + 4.2)/2, None, 1.1),
            # RHC
            ('BPs', max(vars_set['brachial_pressure'][interval])/mmHg2SI, 118, None, 27),
            ('BPd', min(vars_set['brachial_pressure'][interval])/mmHg2SI, 62, None, 3),
            # ('BPk', numpy.mean(vars_set['renal_capillary'][interval]) /mmHg2SI, 20, None, 1, 1/mmHg2SI),
            ('Ppas', numpy.max(vars_set['P_pa'][interval])/mmHg2SI, 38, None, 11),
            ('Ppad', numpy.min(vars_set['P_pa'][interval])/mmHg2SI, 21, None, 7),
            ('Ppv', numpy.mean(vars_set['P_pv'][interval])/mmHg2SI, 17, None, 8),
            ('P_RVs', numpy.max(vars_set['heartComponent.ventricles.P_RV'][interval])/mmHg2SI, 41.6, None, 21),
            ('P_RVd', numpy.min(vars_set['heartComponent.ventricles.P_RV'][interval])/mmHg2SI, 6, None, 8.5),
            ]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], std=o[4], k_p=1), ov))


    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives
