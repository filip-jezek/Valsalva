import scipy.io as skipy
import fun_lib
import re
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function

def primaprim():
    print("ok")

def getValueNames():
    # for the output file header
    # split to array so we can mix it with costs values and %
    s = '(BPs, BPd, BPk, CO, EF, pwv)'
    return s.strip('()').split(', ')

def getCurrentValues(vars_set):
    
    interval = fun_lib.findInterval(25, 30, vars_set['time'])
    BPs = max(vars_set['brachial_pressure'][interval])
    BPd = min(vars_set['brachial_pressure'][interval])
    
    CO = sum(vars_set['CO'][interval]) / len(interval)
    EF = fun_lib.calculateEF(vars_set['V_LV'][interval])
    # using simpler mean
    BPk = sum(vars_set['renal_capillary'][interval]) / len(interval)

    # Van Bortel 2012 siggest using 80 % of carotid to femoral distance
    # distance = vars_set['Systemic1.speedSegmentLength'][1]*0.8
    distance = 0.672*0.8

    pwv = fun_lib.calculatePWV(
        vars_set['time'], 
        vars_set['carotid_pressure'], 
        vars_set['femoral_pressure'], 
        distance)

    return (BPs, BPd, BPk, CO, EF, pwv)


def calculateCosts(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60


    BPs_target = 120*mmHg2SI
    BPd_target = 80*mmHg2SI
    BPk_target = 20*mmHg2SI
    CO_target = 6.3*lpm2SI
    EF_target = 0.6
    pwv_bounds = [2, 10]


    curvals = getCurrentValues(vars_set)

    costs = list(map(fun_lib.cost, 
        curvals[0:-1], 
        (BPs_target, BPd_target, BPk_target, CO_target, EF_target)))
    
    costs.append(fun_lib.penalty(curvals[-1], pwv_bounds[0], pwv_bounds[1]))
   
    return costs