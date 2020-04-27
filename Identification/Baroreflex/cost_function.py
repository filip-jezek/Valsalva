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

    BPs_target = 160*mmHg2SI
    EDV_target = 133*ml2SI
    ESV_target = 27*ml2SI

    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(10, 15, vars_set['time'])
    interval2 = fun_lib.findInterval(10, 15, vars_set['time'])

    # build costs
    ov = [  ('costs', max(vars_set['sum_cost'][interval]), -1, None, 1),
            ('fbr_aor', max(vars_set['systemicMockPressure.baroreflex_system.baroreceptor_aortic.fbr'][interval2]), 35, None, 10),
            ('fbr_car', max(vars_set['systemicMockPressure.baroreflex_system.baroreceptor_carotid.fbr'][interval2]), 35, None, 10),
            ]
    
    # make it a dict?
    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight = o[4]), ov))
    return objectives
