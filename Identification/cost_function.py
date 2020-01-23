import scipy.io as skipy
# import DyMat
# from matplotlib import pyplot as plt

# // calculate the cost function

def WriteCost(vars_set):
    total_cost = repr(calculateCost(vars_set))
    with open('dsout.out', 'w') as file:
        file.write("banik pico\n")
        file.write('f(x) =' + total_cost)
    

def calculateCost(vars_set):

    Pa = vars_set['Systemic1.aortic_arch_C2.port_a.pressure']
    
    interval = findInterval(280, 300, vars_set['time'])
    Pa_avg = sum(Pa[interval]) / len(interval)

    cost = Pa_avg - 100*133.32

    return cost
# // write the outputfiles

def findInterval(t_from, t_to, timeArr):
    return range(findLowestIndex(t_from, timeArr), findLowestIndex(t_to, timeArr))

def findLowestIndex(time, timeArr):
    lst = timeArr.tolist()
    return next((i for i, x in enumerate(lst) if x >= time))
