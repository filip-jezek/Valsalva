import scipy.io as skipy
import fun_lib
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function

def getObjectives(vars_set):

    endTime = vars_set['time'][-1]
    interval = fun_lib.findInterval(endTime - 10, endTime, vars_set['time'])
    BPs = max(vars_set['brachial_pressure'][interval])
    BPd = min(vars_set['brachial_pressure'][interval])
    
    # this shows quite a too low CO, probably due to sampling, sharp peaks and no interpolation
    # CO = sum(vars_set['heartComponent.aorticValve.q_in.q'][interval]) / len(interval)

    CO = sum(vars_set['CO'][interval]) / len(interval)

    objectives = list()

    def buildObjective(name, value, targetVal):
        o = fun_lib.ObjectiveVar(name, value, targetVal)
        objectives.append(o)

    buildObjective('BPs', BPs, 120*133.32)
    buildObjective('BPd', BPd, 80*133.32)
    buildObjective('CO', CO, (6.34/1000/60))

    return objectives