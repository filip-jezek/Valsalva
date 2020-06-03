import numpy
import math
import enum
from scipy.signal import savgol_filter
import numpy
import scipy.signal as ss
import os
import re
from typing import Iterable

# For variance based cost function we need a guess of variance of any target variables
# does not really work for timed variables though
DEFAULT_STD_PERCENT = 10

class CostFunctionType(enum.Enum):
    Ignore = 0
    Quadratic = 1
    Linear = 2
    QuadraticVariance = 3
    DistanceFromZero = 4
    

class ObjectiveVar:
    """
    Provides objects and functions for calculating cost functions

    """

    def __init__(self, name, 
                value=None, 
                targetValue=None, 
                limit=None, 
                weight=1, 
                k_p=1e3, 
                std = None,
                costFunctionType=CostFunctionType.Quadratic):
        self.name = name
        self.targetValue = targetValue
        self.value = value
        self.limit = limit if limit is not None else [-math.inf, math.inf]
        self.weight = weight
        # self.costLimit = -1 # unlimited costs
        self.k_p = k_p # multiplier for out ouf limit values
        self.costFunctionType = costFunctionType
        # standard deviation
        self.std = std

    def __cost_function(self, measured, target):
        # calculate costs. Could go negative or NaN for negative or zero measured values!
        if self.costFunctionType is CostFunctionType.Quadratic:
            return self.weight*(measured - target)**2/(target**2)
        elif self.costFunctionType is CostFunctionType.Linear:
            return self.weight*abs(measured - target)/target
        elif self.costFunctionType is CostFunctionType.QuadraticVariance:
            # variance is squared standard deviation
            return self.weight*(measured - target)**2/(target*self.std)
        elif self.costFunctionType is CostFunctionType.DistanceFromZero:
            return self.weight*self.targetValue
        elif self.costFunctionType is CostFunctionType.Ignore:
            return 0
        else:
            raise NotImplementedError()



    def cost(self):
        """
        Calculates costs per variable as difference from target value and applies steep linear penalty for outside of the bounds values, if defined.
        """
        measured = self.value
        if self.targetValue is None:
            c = 0
        else:
            if self.std is None:
                self.std = self.targetValue*DEFAULT_STD_PERCENT/100

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


def getObjectiveByName(objective_iterable, name):
    """
    Returns objective of given name from the list
    """
    return next((i for i in objective_iterable if i.name == name))

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

def getOddWindow(time, dt):
    win = round(time/dt)
    return int(win) if win%2 == 1 else int(win) + 1

def detrend(sig, window, cutoff = -math.inf):
    sigf = savgol_filter(sig, window, 3)
    return [max(s - sf, cutoff) for s, sf in zip(sig, sigf)]

def getPeaks(sig, dt):
    # get mean and detrend
    win = getOddWindow(2, dt)
    sigDet = detrend(sig, win, 0)
    # and again
    sigDet2 = detrend(sigDet, win, 0)
    # plt.plot(sigDet)
    # plt.plot(sigDet2)
    # plt.show()

    # find peaks - its minimal distance is 0.3 s (is about 200 BPM) and is above the mean of the signal
    peaks, _ = ss.find_peaks(sigDet2, distance= int(0.4/dt))
    return peaks

def getHR(sig, peaks, dt):
    
    heart_rate = [0]*len(sig)
    for i in range(1, len(peaks)):
        # loop from 2nd
        # take range inbetween the means
        rng = slice(peaks[i-1], peaks[i], 1)
        # time between the peaks
        beat_time = (peaks[i] - peaks[i-1])*dt
        heart_rate[rng] = [1/beat_time]*(rng.stop-rng.start)
     
    # before the first peak, lets fill with the first value
    heart_rate[0: peaks[0]] = [heart_rate[peaks[0]]]*peaks[0]

    # fill in the last
    l = len(heart_rate[peaks[-1]:] )
    heart_rate[peaks[-1]:] = [heart_rate[peaks[-1]-1]]*l

    return numpy.array(heart_rate)
    

def getMeanRR(sig, peaks):

    means = [0]*len(sig)
    for i in range(1, len(peaks)):
        # loop from 2nd
        # take range inbetween the means
        rng = slice(peaks[i-1], peaks[i], 1)
        # take mean from two peaks
        means[rng] = [sig[rng].mean()]*(rng.stop-rng.start)
    
    # fill in the begining - value of first peak, long up to first peak indice
    means[0: peaks[0]] = [means[peaks[0]]]*peaks[0]

    # fill in the last
    l = len(means[peaks[-1]:] )
    means[peaks[-1]:] = [means[peaks[-1]-1]]*l

    return numpy.array(means)

def getValsalvaStart(time, thoracic_pressure, threshold = 1330):
    
    return time[findLowestIndex(threshold, thoracic_pressure)]

def getValsalvaEnd(valsalva_start, time, thoracic_pressure, min_length = 5, threshold = 1330):
    # minimal valsalva length is 5
    i_from = findLowestIndex(valsalva_start + min_length, time)
    # get rid of noise for thresholding
    tpf = ss.medfilt(thoracic_pressure, 7)
    
    i = next((i for i, x in numpy.ndenumerate(tpf[i_from:]) if x <= threshold))
    return time[i] + valsalva_start + min_length

def updateObjectivesByValuesFromFile(filename, objectives = None) -> Iterable[ObjectiveVar]:
    """
    If no objectives are provided, get a list of all objectives listed in the file instead
    """
    if objectives is None:
        create = True
        objectives = list()
    else:
        create = False

    with open(filename, 'r') as file:
        lines = file.readlines()

        for line in lines[1:]:
            # first col is name
            vals = line.split(',')
            if create:
                objective = ObjectiveVar(vals[0], costFunctionType=CostFunctionType.Ignore)
                objectives.append(objective)
            else:
                objective = next((o for o in objectives if o.name == vals[0]), None)

            if objective is not None:
                objective.targetValue = float(vals[1])
                objective.std = float(vals[2])
    
    return objectives

def getRunNumber():
    """ Gets GenOpt run number using the name of the current working directory
    """
    cur_dirname = os.path.basename(os.getcwd())
    run_match = re.match(r'[\w-]*-(\d+)$', cur_dirname)

    if run_match is not None:
        return int(run_match[1])
    else:
        return 0

def getSafeLogDir(unsafeDir):
    """ Try provided unsafeDir and falls back to current dir otherwise
    """

    if not os.path.isdir(unsafeDir):
        return ''
    else:
        return unsafeDir + '\\'

