import numpy
import math

class ObjectiveVar:

    def __init__(self, name, value = None, targetValue = None, limit = None, weight = 1, k_p = 1e3):
        self.name = name
        self.targetValue = targetValue
        self.value = value
        self.limit = limit if limit is not None else [-math.inf, math.inf]
        self.weight = weight
        # self.costLimit = -1 # unlimited costs
        self.k_p = k_p # multiplier for out ouf limit values

    def __cost_function(self, measured, target):
        # calculate costs. Could go negative or NaN for negative or zero measured values!
        return self.weight*(measured - target)**2/(measured*target)

    def cost(self):
        """
        Calculates costs per variable as difference from target value and applies steep linear penalty for outside of the bounds values, if defined.
        """
        measured = self.value
        if self.targetValue is None:
            c = 0
        else:
            c = self.__cost_function(measured, self.targetValue)

        min_val = self.limit[0]
        max_val = self.limit[1]

        if measured < min_val:
            return c + self.__cost_function(measured, min_val)*self.k_p
        elif measured > max_val:
            return c + self.__cost_function(measured, max_val)*self.k_p
        else:
            return c
    
    def inLimit(self):
        if self.value < self.limit[0] or self.value > self.limit[1]:
            return False
        else:
            return True


    
# Provides objects and functions for calculating cost functions

def findInterval(t_from, t_to, timeArr):
    return range(findLowestIndex(t_from, timeArr), findLowestIndex(t_to, timeArr))

def findLowestIndex(time, timeArr):
    lst = timeArr.tolist()
    return next((i for i, x in enumerate(lst) if x >= time))


def avg(var_set, interval):
    return sum(var_set[interval]) /len(interval)

def calculatePWV(timeArr, signal, delayedSignal, distance):
    """
    Calculates pulse wave propagation velocity based on bottom base shift of pressure signals.
    In other words, it compares  position of signals minimum

    This is the simpler approach, however it does not work in all cases.
    """    
    # find last three seconds - should do for as low as 30 bpm
    # we better leave some time to find min of the other signal
    i_from = findLowestIndex(timeArr[-1] - 3, timeArr)
    i_to = findLowestIndex(timeArr[-1] - 1, timeArr)

    # signal window
    sig_win = signal[i_from:i_to].tolist()
    i_first_min = sig_win.index(min(sig_win)) + i_from
    t_first_min = timeArr[i_first_min]

    # 2 m/s is an absolute lazy crazy slow velocity to cut the second signal window with
    max_t = distance / 2
    i_second_window = findLowestIndex(t_first_min + max_t, timeArr)
    # second window is shifted from the first one and substantially shorter
    delayedSignal_window = delayedSignal[i_first_min:i_second_window].tolist()
    i_second_min = i_first_min + delayedSignal_window.index(min(delayedSignal_window))
    t_second_min = timeArr[i_second_min]

    timediff = t_second_min - t_first_min
    velocity = distance/timediff
    return velocity

def calculatePWV2(timeArr, signal, delayedSignal, distance):
    """
    More robust way to calculate PWV takes position of maximal signal rising steepness.
    In other words, takes maxs of difference

    """

    diff_signal = numpy.diff(signal)
    diff_delayedSignal = numpy.diff(delayedSignal)
    
    # call with signal derivatives - maxing the steepness
    # we have to use negative though, as calculatePWV is taking mins instead of max
    return calculatePWV(timeArr, numpy.negative(diff_signal), numpy.negative(diff_delayedSignal), distance)

def calculateEF(volumes):
    esv = min(volumes)
    edv = max(volumes)
    return (edv - esv)/edv

