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

    BPs_target = 157
    BPd_target = 87
    BPm_target = 114
    BPk_target = 20
    

    interval = fun_lib.findInterval(25, 30, vars_set['time'])

    BPs = max(vars_set['brachial_pressure'][interval])
    BPd = min(vars_set['brachial_pressure'][interval])

    # using true mean
    BPm = sum(vars_set['brachial_pressure_mean'][interval]) / len(interval)

    # using simpler mean
    BPk = sum(vars_set['renal_capillary'][interval]) / len(interval)

    print('BPs %s mmHg (%s)' % tuple(map(lambda x: str(round(x)), (BPs/mmHg2Pa, BPs_target))))
    print('BPd %s mmHg (%s)' % tuple(map(round, (BPd/mmHg2Pa, BPd_target))))
    print('SV %s  ml (%s)' % tuple(map(round, (BPm/mmHg2Pa, BPm_target))))

    cost = lambda measured, target: (measured / target - 1)**2
    costs = list(map(cost, (BPs, BPd, BPm, BPk), (BPs_target*mmHg2Pa, BPd_target*mmHg2Pa, BPm_target*mmHg2Pa, BPk_target*mmHg2Pa)))
   

    # costs = [(BPs / (BPs_target*mmHg2Pa) - 1)**2,  (BPd / (80*133.32) - 1)**2, (CO / (6.34/1000/60) - 1)**2]

    return costs