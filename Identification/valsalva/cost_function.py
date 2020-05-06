import scipy.io as skipy
import fun_lib
from statistics import mean
from scipy.signal import savgol_filter
from matplotlib import pyplot as plt
import scipy.signal as ss
import math
# import DyMat
# from matplotlib import pyplot as plt

def getOddWindow(time, dt):
    win = round(time/dt)
    return int(win) if win%2 == 1 else int(win) + 1


def getMeanSavGol(sig, dt):
    # smooth the data first to have the mean
    # we assume equidistant grid
    # dt = vars_set['time'][2] - vars_set['time'][1]
    # dt = time[2] - time[1]
    sig_mean = savgol_filter(sig, getOddWindow(3, dt), 3) # window size of 3 following avg heartbeats = 3s and polynomial order of 3

    return sig_mean

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

def getMeanRR(sig, dt):
    
    peaks = getPeaks(sig, dt)
    means = [0]*len(peaks)
    for i in range(2, len(peaks)):
        # loop from 2nd
        # take range inbetween the means
        rng = slice(peaks[i-1], peaks[i], 1)
        # take mean from two peaks
        means[rng] = [sig[rng].mean()]*(rng.stop-rng.start)
    return means

def getObjectives(vars_set):

    # Pa = vars_set['Systemic#1.aortic_arch_C2.port_a.pressure']
    # Pa = vars_set['Pa']
    # interval = findInterval(380, 400, vars_set['time'])
    # Pa_avg = sum(Pa[interval]) / len(interval)

    # HR is system input (change in phi) from 70 to 89
    mmHg2SI = 133.32
    ml2SI = 1e-6
    lpm2SI = 1e-3/60

    BP = vars_set['brachial_pressure']
    dt = vars_set['time'][2] - vars_set['time'][1]
    assert dt > 0, "The simulation must not store values at events."
    bp_mean = getMeanRR(BP, dt)
    time = vars_set['time']

    bp_mean1 = getMeanSavGol(BP, dt)
    bp_mean2 = getMeanRR(BP, dt)
    # plt.plot(time, BP)
    # plt.plot(time, bp_mean1)
    # plt.plot(time[0:len(bp_mean2)], bp_mean2)
    # sigf = ss.medfilt(BP, int(1/dt)+1)
    # plt.plot(time[0:len(sigf)], sigf)

    # plt.legend(['bp', 'mean1', 'mean2', 'medfilt'])
    
    # plt.show()


    # find valsalva start and end
    valsalva_start = 20
    # TODO find automatically
    valsalva_end = 35

    # baseline is 5s before Valsalva
    baseline = bp_mean[fun_lib.findInterval(valsalva_start-5, valsalva_start)].mean()

    t = vars_set['time'][-1]
    interval = fun_lib.findInterval(valsalva_start, valsalva_end, vars_set['time'])
    interval2 = fun_lib.findInterval(10, 15, vars_set['time'])


    # systolic peak at start of valsalva
    interval_phase1 = fun_lib.findInterval(valsalva_start, valsalva_start + 5, vars_set['time'])
    ph1_peak = max(BP[interval])
    
    # min mean pressure during valsalva

    interval_valsalva = fun_lib.findInterval(valsalva_start, valsalva_end, vars_set['time'])
    valsalva_mean_min = min(bp_mean1[interval_valsalva])
    # peak minimal after valsalva release
    # peak maximal after valsalva
    # mean recovery


    # build costs
    ov = [  ('costs', max(vars_set['sum_cost'][interval]), -1, None, 1),
            ('fbr_aor', max(vars_set['systemicMockPressure.baroreflex_system.baroreceptor_aortic.fbr'][interval2]), 35, None, 10),
            ('fbr_car', max(vars_set['systemicMockPressure.baroreflex_system.baroreceptor_carotid.fbr'][interval2]), 35, None, 10),
            ]
    
    # make it a dict?
    objectives=list(map(lambda o: fun_lib.ObjectiveVar(o[0], value = o[1], targetValue = o[2], limit=o[3], weight = o[4]), ov))
    return objectives
