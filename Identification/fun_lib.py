import numpy
# Provides functions for calculating cost functions

def findInterval(t_from, t_to, timeArr):
    return range(findLowestIndex(t_from, timeArr), findLowestIndex(t_to, timeArr))

def findLowestIndex(time, timeArr):
    lst = timeArr.tolist()
    return next((i for i, x in enumerate(lst) if x >= time))

def cost(measured, target):
    # calculate costs. Could go negative or NaN for negative or zero measured values!
    return (measured - target)**2/(measured*target)

def penalty(measured, min_val, max_val, k_p = 1e3):
    # apply steep linear penalty for outside of the bounds values
    if measured < min_val:
        return cost(measured, min_val)*k_p
    elif measured > max_val:
        return cost(measured, max_val)*k_p
    else:
        return 0

def calculatePWV(timeArr, signal, delayedSignal, distance):
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

def calculateEF(volumes):
    esv = min(volumes)
    edv = max(volumes)
    return (edv - esv)/edv

