import scipy.io as skipy
import matplotlib.pyplot as plt
import fun_lib
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function

def plotObjectives(vars_set, interval, objectives):
    if '__plot_axes' in vars_set:
        ax = vars_set['__plot_axes']
    else:
        fig = plt.figure()
        ax = fig.subplots()

    ax.plot(vars_set['time'], vars_set['brachial_pressure']/133.32, label='Brachial pressure mmHg')
    ax.plot(vars_set['time'], vars_set['CO']*1000*60, label='CO l/min')

    # bounds
    pack = (objectives, vars_set, ax, interval)
    fun_lib.plotObjectiveTarget(pack,'BPs', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'BPd', 1/133.32)
    fun_lib.plotObjectiveTarget(pack,'CO', 1000*60, fmt='r')

    total_costs = fun_lib.countTotalWeightedCost(objectives)
    ax.set_title('Baseline costs %.6f' % total_costs)


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

    # to have comparable cost function values one must have the stds ready
    map(fun_lib.unifyCostFunc, objectives)

    if '__draw_plots' in vars_set and vars_set['__draw_plots']:
        plotObjectives(vars_set, interval, objectives)

    return objectives