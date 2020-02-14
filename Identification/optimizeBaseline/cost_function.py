import scipy.io as skipy
import fun_lib
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function

def primaprim():
    print("ok")

def calculateCosts(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # cost = (Pa_avg - 100*133.32)**2

    interval = fun_lib.findInterval(50, 60, vars_set['time'])
    BPs = max(vars_set['Systemic1.brachial_L82.u_C'][interval])
    BPd = min(vars_set['Systemic1.brachial_L82.u_C'][interval])
    
    # this shows quite a too low CO, probably due to sampling, sharp peaks and no interpolation
    # CO = sum(vars_set['heartComponent.aorticValve.q_in.q'][interval]) / len(interval)

    CO = sum(vars_set['heartComponent.CO_LV2'][interval]) / len(interval)

    costs = [(BPs / (120*133.32) - 1)**2,  (BPd / (80*133.32) - 1)**2, (CO / (6.34/1000/60) - 1)**2]

    return costs