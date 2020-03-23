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
    s = '(EF)'
    return s.strip('()').split(', ')

def getCurrentValues(vars_set):
    
    interval = fun_lib.findInterval(25, 30, vars_set['time'])
    EDV = max(vars_set['V_LV'][interval])
    ESV = min(vars_set['V_LV'][interval])
    
    SV = EDV - ESV
    EF = fun_lib.calculateEF(vars_set['V_LV'][interval])

    return [EF]


def calculateCosts(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60


    EF_target = 0.6
    SV_target = 90*ml2SI


    curvals = getCurrentValues(vars_set)

    # costs = list(map(fun_lib.cost, 
    #     curvals, 
    #     [EF_target]))
    
    costs = [curvals[0] - EF_target]
       
    return costs