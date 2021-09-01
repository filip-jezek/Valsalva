import numpy
import math
import enum
from scipy.signal import savgol_filter
import numpy
import scipy.signal as ss
import os
import re
import matplotlib.pyplot as plt
from typing import Iterable
import importlib.util
import inspect

# For variance based cost function we need a guess of variance of any target variables
# does not really work for timed variables though
DEFAULT_STD_PERCENT = 10

class SimulationResultIncompleteError(ValueError):
    def __init__(self, length:float, minimalLength:float, owner:str, penalty):
        m = 'Incomplete simulation output, the simulation probably crashed at %.2f s (%.1f s required). At %s' % (length, minimalLength, owner)
        super().__init__(m)
        self.penalty = penalty

class CostFunctionType(enum.Enum):
    Ignore = 0
    Quadratic = 1
    QuadraticVariance = 2
    Linear = 3
    LinearVariance = 4
    DistanceFromZero = 5
    

class ObjectiveVar:
    """
    Provides objects and functions for calculating cost functions

    """

    def __init__(self, name, 
                value=None, 
                targetValue=None, 
                limit=None, 
                weight=None,
                tolerance = None,
                k_p=1e3, 
                std = None,
                base = 1,
                costFunctionType=CostFunctionType.Quadratic):
        self.name = name
        self.targetValue = targetValue
        self.value = value
        self.limit = limit if limit is not None else [-math.inf, math.inf]

        if weight is None:
            if tolerance is not None and targetValue is None and limit is not None:
                # calculate the weight from limits
                self.weight= 0 if tolerance == 0 else max(limit)/tolerance 
            elif tolerance is not None and targetValue is not None:
                # calculate the weigth from tolerance and target value
                self.weight = 0 if tolerance == 0 else targetValue/tolerance
            else:
                # cant calculate the weight
                self.weight = 1
        else:
            self.weight = weight

        # self.costLimit = -1 # unlimited costs
        self.k_p = k_p # penalty multiplier for out ouf limit values
        self.costFunctionType = costFunctionType
        # standard deviation
        self.std = std
        self.base_cost = None
        self.base = base

    def __cost_function(self, measured, target):
        # calculate costs. Could go negative or NaN for negative or zero measured values!
        if self.costFunctionType is CostFunctionType.Quadratic:
            if target is None:
                return 0
            elif target == 0:
                return self.weight*abs(measured)
            else:
                return self.weight*(measured - target)**2/(target**2)

        elif self.costFunctionType is CostFunctionType.QuadraticVariance:
            if target is None:
                return 0            
            elif target == 0:
                return self.weight*abs(measured)
            if self.std is None:
                self.std = abs(target*DEFAULT_STD_PERCENT/100)
            
            # variance is squared standard deviation
            return self.weight*(measured - target)**2/(target*self.std)

        elif self.costFunctionType is CostFunctionType.Linear:
            if target is None:
                return 0
            elif target == 0:
                return self.weight*abs(measured)
            
            return self.weight*abs(measured - target)/abs(target)

        elif self.costFunctionType is CostFunctionType.LinearVariance:
            if target is None:
                return 0
            elif target == 0:
                return self.weight*abs(measured)
            if self.std is None:
                self.std = abs(target*DEFAULT_STD_PERCENT/100)
            # variance is squared standard deviation
            return self.weight*abs(measured - target)/(self.std)

        elif self.costFunctionType is CostFunctionType.DistanceFromZero:
            return self.weight*abs(measured)

        elif self.costFunctionType is CostFunctionType.Ignore:
            return 0
        else:
            raise NotImplementedError()



    def cost(self):
        """
        Calculates costs per variable as difference from target value and applies steep linear penalty for outside of the bounds values, if defined.
        """
        measured = self.value
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

    def target_log(self, normalizeToBase = False):
        if normalizeToBase:
            base = self.base
        else:
            base = 1.0

        if self.costFunctionType is CostFunctionType.DistanceFromZero:
            target = "%d" % 0
        elif self.targetValue is not None:
            target = "%.3g" % (self.targetValue*base)
        elif self.limit is not None:
            if self.value < self.limit[1]:
                # lower than upper bound - print the lower one
                target = '>%0.3g' % (self.limit[0]*base)
            else: # self.value >
                target = '<%0.3g' % (self.limit[1]*base)
                # target = ' in limit' if self.inLimit() else 'out limit'
        else:
            target = ' ??wat?? '
        
        return target.ljust(6)

    def __repr__(self):
        return '%s,%.3g,%s' % (self.name, self.value, self.target_log())

    def cost_to_base_comparisson(self) -> float:
        c = self.cost()
        if self.base_cost is None:
            return ''
        elif self.base_cost == 0 and c == 0:
            return ''
        elif self.base_cost == 0 and c > 0:
            return '+Inf'
        elif self.base_cost == 0 and c < 0:
            # improbable :)
            return '-Inf'
        else:
            return '%+0.0g' % ((self.cost() - self.base_cost)/self.base_cost*100)


def getObjectiveByName(objective_iterable, name) -> ObjectiveVar:
    """
    Returns objective of given name from the list
    """
    try:
        return next((i for i in objective_iterable if i.name == name))
    except StopIteration:
        return None

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

def calculateQdot_mv(v_lv, q_mv, time, interval):

    # find the v_lv min in first half of the interval
    half_interval = range(interval[0], int((interval[-1] - interval[0])/2 + interval[0]))
    # end systolic time index
    t_ES_i = numpy.argmin(v_lv[half_interval]) + half_interval[0]
    
    # security countermeasure for weird flor profiles
    buffer = next(i for  i, v in enumerate(q_mv[t_ES_i:]) if v > 0)
    if buffer is not None and buffer > 0 and buffer < len(half_interval):
        t_ES_i = t_ES_i + buffer

    ESV = v_lv[t_ES_i]

    # difference in the flow - positive for increasing flow rate, negative for decreasing flow
    dq_mv = numpy.diff(q_mv)

    # find the index of just after first mitral flow peak - i.e. the first occurence of slowing the flow down.
    # ignore zero, as this might mean an event
    t_q_mv_peak1_i = next(i for  i, v in enumerate(dq_mv[t_ES_i:]) if v < 0) + t_ES_i

    # from t_q_mv_peak1 find first increase in rate of flow - thats the beginning of the atrial kick
    t_q_mv_saddle_i = next(i for  i, v in enumerate(dq_mv[t_q_mv_peak1_i:]) if v > 0) + t_q_mv_peak1_i

    # and thats the pre-atrial-kick volume we are interested in
    V_LV_passive = v_lv[t_q_mv_saddle_i]

    # find EDV time index as first zero or negative mitral flow AFTER the first peak
    t_ED_i = next(i for  i, v in enumerate(q_mv[t_q_mv_peak1_i:]) if v <= 0) + t_q_mv_peak1_i
    EDV = v_lv[t_ED_i]

    # atrial kick fraction as a fractions of volumes
    SV = (EDV - ESV)
    if SV == 0:
        atrial_kick_fraction = 0
    else:
        atrial_kick_fraction = (EDV - V_LV_passive)/SV

    q_mv_saddle = q_mv[t_q_mv_saddle_i]


    print('Atrial kick %d%% from %1.2e m3/s, timings: ESV %dml at %.2fs, passive %dml at %.2fs, EDV %dml at %.2fs' % 
        (round(atrial_kick_fraction*100), 
        q_mv_saddle,
        round(v_lv[t_ES_i]*1e6), 
        time[t_ES_i], 
        round(v_lv[t_q_mv_saddle_i]*1e6), 
        time[t_q_mv_saddle_i], 
        round(v_lv[t_ED_i]*1e6), 
        time[t_ED_i]))

    return (atrial_kick_fraction, q_mv_saddle)




def calculateQ_MV(q, time, interval):
    """ Returns tuple of passive and atrial kick peaks mitral flow heart filling rate

    q : flow signal
    time : time set
    interval : guess range interval in which to search. May overflow during search, so leave at least one beat buffer.

    returns
    """
    # find the first peak in the interval
    max_index = numpy.argmax(q[interval]) + interval[0]

    # find an one beat interval - lloking for first zero or negative flow (closed valve)
    min_index = next(i for  i, v in enumerate(q[max_index:]) if v <= 0) + max_index

    # get the peak in reverse direction
    m = 0
    for v in reversed(q[max_index:min_index]):
        # find a maximum
        m = max(m, v)
        # detect a decrease by 10% from max and stop right there, before hitting the main peak
        if v <= 0.9*m:
            break
    
    return (q[max_index], m)

def calculateEA(q_mv,cc, time, interval, A_length = 0.2):
    """ Returns tuple of passive and atrial kick peaks mitral flow heart filling rate
    assumes the passive peak is larger

    q : flow signal
    cc : cardiac cycle phase to get the systole from. Cardiac cycle starts with atrial contraction.
    time : time set
    interval : guess range interval in which to search. May overflow during search, so leave at least one beat buffer.

    returns
    """
    
    # find begging of the cardiac cycle - that is the beggingin of atrial systole
    A_start_i = numpy.argmax(cc[interval]) + interval[0]
    
    # find index of time + 0.2
    A_end_i = numpy.argmax(time > time[A_start_i] + A_length)
    A_i = numpy.argmax(q_mv[A_start_i:A_end_i]) + A_start_i
    A = q_mv[A_i]
    # find the saddle - its between start of the cycle and maximal A flow
    q_sad = numpy.min(q_mv[A_start_i:A_i + 1])

    # find the max peak in the interval - that is the E peak
    # E_max_i = numpy.argmax(q_mv[interval]) + interval[0]
    E = numpy.max(q_mv[interval])



    # ratio
    if A != 0:
        EA = E/A
    else:
        EA = 0
    
    if q_sad != 0:
        As = A/q_sad
    else:
        As = 0

    return (EA, As)

def calculateAlternativeCO(vars_set, interval):
    """ Workaround for when the CO is not available"""

    if 'CO' not in vars_set or (vars_set['CO'][interval]).all() == 0:
        co = (numpy.max(vars_set['V_LV'][interval]) - numpy.min(vars_set['V_LV'][interval]))*numpy.mean(vars_set['HR'][interval])
    else:
        co = numpy.mean(vars_set['CO'][interval])

    return co

    

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

def getPPulseRR(sig, peaks):
    """ Returns pulse pressure inbetween the peaks """

    pp = [0]*len(sig)
    for i in range(1, len(peaks)):
        # loop from 2nd
        # take range inbetween the means
        rng = slice(peaks[i-1], peaks[i], 1)
        # take mean from two peaks
        pp[rng] = [sig[rng].max() - sig[rng].min()]*(rng.stop-rng.start)
    
    # fill in the begining - value of first peak, long up to first peak indice
    pp[0: peaks[0]] = [pp[peaks[0]]]*peaks[0]

    # fill in the last
    l = len(pp[peaks[-1]:] )
    pp[peaks[-1]:] = [pp[peaks[-1]-1]]*l

    return numpy.array(pp)

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
        print('Getting objectives from %s ...' % filename)
    else:
        create = False
        print('Updating objectives from %s ...' % filename)

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

def getRunNumber() -> str:
    """ Gets GenOpt run number using the name of the current working directory
    """
    cur_dirname = os.path.basename(os.getcwd())
    run_match = re.match(r'[\w-]*-(\d+)$', cur_dirname)

    if run_match is not None:
        return int(run_match[1])
    else:
        return 0

def getSafeLogDir(unsafeDir, safe_rollback = '..\\'):
    """ Try provided unsafeDir and falls back to parent dir otherwise
    """
    if len(unsafeDir) == 0:
        # its current dir, thats fine
        return ''
    
    unsafeDir_path = unsafeDir.rstrip('\\') + '\\'

    if not os.path.isdir(unsafeDir_path):
        return safe_rollback
    else:
        return unsafeDir_path

def unifyCostFunc(o:ObjectiveVar):
    """ Changes the objective cost function type to variance
    """
    if o.costFunctionType == CostFunctionType.Linear:
        o.costFunctionType = CostFunctionType.LinearVariance
    elif o.costFunctionType == CostFunctionType.Quadratic:
        o.costFunctionType = CostFunctionType.QuadraticVariance

def countTotalWeightedCost(objectives : Iterable[ObjectiveVar]) -> float:
    """ Returns total objective costs, weighted by number of objectives - effectively returns mean of costs"""
    active_obj = sum(1 for o in objectives if o.costFunctionType is not CostFunctionType.Ignore)
    
    total_cost = sum(o.cost() for o in objectives)
    # weighing cost by numbr of objectives
    return total_cost / active_obj

def countTotalSumCost(objectives : Iterable[ObjectiveVar]) -> float:
    "Returns simple sum of all objectives costs"
    return sum(o.cost() for o in objectives)    

def getAxes(vars_set : dict) -> plt.axes:
    """ Gets axes from vars_set if defined, creates empty subplots figure otherwise 
    """

    if '__plot_axes' in vars_set:
        return vars_set['__plot_axes']
    else:
        fig = plt.figure()
        return fig.subplots()

def plotObjectiveTarget(pack:tuple, objective_name:str, unitFactor:float, fmt = 'k', verticalalignment = 'bottom', showVal = True):
    """ Plots the objective target with label
    pack = (objectives:list, ax:plt.axes, interval:range)
    """
    (objectives, time, ax, interval) = pack

    try:
        # get the bounds from the target value
        objective = next((o for o in objectives if o.name == objective_name))
    except StopIteration:
        # defensive programming - not all objectives are actrive rn
        return

    if objective.targetValue is None:
        # defensive programming: it might be none
        return
    targetVal = objective.targetValue*unitFactor
    val = objective.value*unitFactor
    ax.plot((time[interval[0]], time[interval[-1]]), [targetVal]*2, fmt)
    if showVal:
        s = '%s = %2f (%.6f)' % (objective_name, val, objective.cost())
    else:
        s = '%s %.6f' % (objective_name, objective.cost())
    ax.text(time[interval[-1]], targetVal, s, 
            horizontalalignment='right', 
            verticalalignment=verticalalignment)

def plotObjectiveLimit(pack:tuple, objective_name:str, unitFactor:float, limit:str, fmt = 'k', verticalalignment = 'bottom', showVal = True):
    """ Plots the objective limit with label
    pack = (objectives:list, time:iterable, ax:plt.axes, interval:range)
    limit = 'lower' or 'upper'
    """
    (objectives, time, ax, interval) = pack
    SV_min_objective = getObjectiveByName(objectives, objective_name)
    limit_val = SV_min_objective.limit[0] if limit == 'lower' else SV_min_objective.limit[1]
    val = SV_min_objective.value*unitFactor
    ax.plot([time[interval[0]], time[interval[-1]]], [limit_val]*2, 'r')
    if showVal:
        s= '%s = %3f, %s lim at %3f (%.4f)' % (objective_name, val, limit, limit_val, SV_min_objective.cost())
    else:
        s = '%s %s lim %.4f' % (objective_name, limit, SV_min_objective.cost())
    ax.text(time[interval[-1]], limit_val, s, horizontalalignment='right', 
            verticalalignment=verticalalignment, fontsize = 8, color='red')

def importCostFunction(dir = '..\\'):
    """ imports cost_function.py from parent directory
    """

    spec = importlib.util.spec_from_file_location(
        'cost_function', dir + 'cost_function.py')
    cf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cf)
    return cf

def checkSimulationLength(simulationTime, minimalSimulationTime, penalty = 1e3):

    if simulationTime < minimalSimulationTime:
        try:
            # this must be errorproof
            frame = inspect.stack()[1]
            filename = frame[0].f_code.co_filename
        except:
            filename = 'Unknwon'

        raise SimulationResultIncompleteError(simulationTime, minimalSimulationTime, filename, penalty)
        
def writeToFile(filename, time:Iterable, signal:Iterable):
    with open(filename, 'w') as file:
        file.write('Time, var\n')
        for t, s in zip(time, signal):
            file.write('%.2f, %.6e\n' % (t, s))
    
    print("Written to %s, mate" % filename)

def logObjectives(objectivesLog_path, objectives, sortBy = 'cost', compare_to_log_path = None, printToStdin=False):
    if compare_to_log_path is not None:
        with open(compare_to_log_path) as file:
            for line in file.readlines():
                # find in objectives
                o = getObjectiveByName(objectives, line.split(',')[0].strip(' '))
                if o is not None:
                    o.base_cost = float(line.split(',')[3])

    if objectivesLog_path is not None:
        if '%' in objectivesLog_path:
            # update with Run number, of supported
            objectivesLog_path = objectivesLog_path % getRunNumber()

        with open(objectivesLog_path, 'w') as file:
            total_cost = countTotalSumCost(objectives)
            if total_cost == 0:
                # fix for division by zero
                total_cost = 1

            max_name_len = max(len(o.name) for o in objectives)

            # tab spaces fro longest header, otherwise for 9 places
            header_names = 'Name'.ljust(max_name_len)


            header = '%s, value    , target   , cost     ,%% tot, [val, target, weight], total = %.6e\n' % (header_names, total_cost)
            file.write(header)
            if printToStdin:
                print(header)

            def sortedObjectives() -> Iterable[ObjectiveVar]:
                if sortBy == 'costs':
                    return sorted(objectives, key = lambda o: o.cost(), reverse=True)
                elif sortBy == 'id':
                    return objectives 
                elif sortBy == 'name':
                    return sorted(objectives, key = lambda o: o.name, reverse=False)
            
            for o in sortedObjectives():
                s = '%s,%.3e ,%s ,%.3e , %02.0f  , %g, %s, %s, %s\n' % (o.name.ljust(max_name_len), o.value, o.target_log(), o.cost(), round(o.cost()/total_cost*100), round(o.value*o.base, 4), o.target_log(normalizeToBase=True), o.cost_to_base_comparisson(), round(o.weight, ndigits=3))
                file.write(s)
                if printToStdin:
                    print(s)
