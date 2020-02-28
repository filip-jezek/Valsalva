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


    # HR is system input (change in phi) from 70 to 89
    mmHg2Pa = 133.32
    ml2m3 = 1e-6
    BPs_target = 113
    BPd_target = 86
    SV_target = 90*0.6 # 60 % of base
    

    interval = fun_lib.findInterval(25, 30, vars_set['time'])

    BPs = max(vars_set['Systemic1.brachial_L82_HeartLevel.u_C'][interval])
    BPd = min(vars_set['Systemic1.brachial_L82_HeartLevel.u_C'][interval])
    
    # this shows quite a too low CO, probably due to sampling, sharp peaks and no interpolation
    # CO = sum(vars_set['heartComponent.aorticValve.q_in.q'][interval]) / len(interval)
    SV = sum(vars_set['heartComponent.SV_LV'][interval]) / len(interval)

    print('BPs %s mmHg (%s)' % tuple(map(lambda x: str(round(x)), (BPs/mmHg2Pa, BPs_target))))
    print('BPd %s mmHg (%s)' % tuple(map(round, (BPd/mmHg2Pa, BPd_target))))
    print('SV %s  ml (%s)' % tuple(map(round, (SV/ml2m3, SV_target))))

    cost = lambda measured, target: (measured / target - 1)**2
    costs = list(map(cost, (BPs, BPd, SV), (BPs_target*mmHg2Pa, BPd_target*mmHg2Pa, SV_target*ml2m3)))
   

    # costs = [(BPs / (BPs_target*mmHg2Pa) - 1)**2,  (BPd / (80*133.32) - 1)**2, (CO / (6.34/1000/60) - 1)**2]

    return costs