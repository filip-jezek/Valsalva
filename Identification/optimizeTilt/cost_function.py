import scipy.io as skipy
import fun_lib
from statistics import mean
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

    BPs_target = 113*mmHg2SI
    BPd_target = 90*mmHg2SI
    SV_target = 63*ml2SI
    phi_target = 0.3739
    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(t-5, t, vars_set['time'])

    # build costs
    ov = [  ('BPs', max(vars_set['brachial_pressure'][interval]), BPs_target, None, 1),
            ('BPd', min(vars_set['brachial_pressure'][interval]), BPd_target, None, 0.1),
            ('SV', fun_lib.avg(vars_set['SV'], interval), SV_target, None, 0.1),
            ('phi', fun_lib.avg(vars_set['phi_baro'], interval), phi_target, None, 1)]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight = o[4]), ov))

    return objectives
