import scipy.io as skipy
import fun_lib
import re
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function


def getObjectives(vars_set):

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
    pwv_bounds = [3.3, 10]

    interval = fun_lib.findInterval(25, 30, vars_set['time'])

    # Van Bortel 2012 siggest using 80 % of carotid to femoral distance
    # distance = vars_set['Systemic1.speedSegmentLength'][1]*0.8
    distance = 0.672*0.8
    pwv = fun_lib.calculatePWV2(
        vars_set['time'],
        vars_set['carotid_pressure'],
        vars_set['femoral_pressure'],
        distance)

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval]), BPs_target, None),
            ('BPd', min(vars_set['brachial_pressure'][interval]), BPd_target, None),
            ('CO', sum(vars_set['CO'][interval]) / len(interval), CO_target, None),
            ('EF', fun_lib.calculateEF(vars_set['V_LV'][interval]), EF_target, None),
            ('BPk', sum(vars_set['renal_capillary'][interval]) / len(interval), BPk_target, None),
            ('PWV', pwv, None, pwv_bounds)            ]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3]), ov))


    return objectives
