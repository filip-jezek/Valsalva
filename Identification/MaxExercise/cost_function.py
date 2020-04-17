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

    EDV_target = 133*ml2SI
    ESV_target = 27*ml2SI

    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(t-5, t, vars_set['time'])

    # build costs
    ov = [  ('EDV', max(vars_set['V_LV'][interval]), EDV_target, None, 1),
            ('ESV', min(vars_set['V_LV'][interval]), ESV_target, None, 1),
            ('Ts', max(vars_set['TS'][interval]), 0.15, None, 1)]

    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight = o[4]), ov))

    return objectives
