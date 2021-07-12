import DyMat
import numpy as np
import time as timer
import os

# experiment_type = 'CHF_baro'
# exp_range = range(100, 0,-5)
# imp_coeff = 1

experiment_type = 'CHF_baro'
exp_range = range(5, 100,5)
imp_coeff = 1

# experiment_type = 'HFpEFDilated_baro'
# exp_range = range(100, 190,10)
# imp_coeff = 1

# experiment_type = 'HFpEF_baro'
# exp_range = range(4700, 5000,100)
# imp_coeff = 1/100

# experiment_type = 'HFpEF_baroInitVol1000'
# exp_range = range(2600, 2700,100)
# imp_coeff = 1/100

filePattern = R"CVS_renalRegulation_%s_%0d.mat"
# filePattern = R"D:\Data\CVS_renalRegulation_CHF_%0d.mat"
# filename = R"c:\home\UMICH\Valsalva\Results2\CardiovascularSystem.mat"


mmHg2SI = 133.322
L2SI = 1e-3
LpD2SI = 1/(1000*60*60*24)

bpm_base = 99*mmHg2SI
startTime = 600


s = '{f_LV}, {vol}, {time}, {eGFR}, {p_int}, {HR}, {CO}, {CVP}, {PCWP}, {EF}, {EDV}, {BPM}, {Qlymph}, {Vint_r},{Vint_excess}, {ESP}, {EDP}, {comp}\n'
result_set = []

outputFile = 'VolumeLoading_%s.csv' % experiment_type
with open(outputFile, 'a') as file:
# with open('CHF_VolumeLoading.csv', 'w') as file:    
    file.write(s.format(f_LV = 'f_LV', vol = 'vol [L]', time = 'time [s]', eGFR = 'eGFR', p_int = 'p_int', HR = 'HR', CO = 'CO [L/min]', EF = 'EF', EDV = 'EDV ml', CVP = 'CVP', PCWP = 'PCWP (mmHg)', BPM = 'maxBP', Qlymph = 'Qlymph [L/day]', Vint_r= 'Relative change to interstitial volume', Vint_excess = 'Excess interstitial volume [L]', ESP = 'ESP', EDP = 'EDP', comp = 'compensated'))

    for i in exp_range:
        filename = filePattern % (experiment_type, i)
        # filename = R"c:\home\UMICH\Valsalva\Results2\CardiovascularSystem.mat"

        try:
            print("Loading file %s..." % filename, end = '')
            tic = timer.time()
            datafile = DyMat.DyMatFile(filename)
            print("Ok [%ds]. Processing..." % (timer.time() - tic), end = '')
        except FileNotFoundError:
            print(" File %0d  not found !" % i)
            continue
        except:
            print(" File %0d  failed, show must go on" % i)

            # os.system("c:\\Program Files\\Dymola 2021\\bin\\dsres2dsf %s %s" % (filename, filename.replace('.mat', '.sdf'))

            continue


        time = datafile.abscissa(2)[0]
        i_dt = int(6/(time[100]/100)) # index of mean interval of 6s , measured over range of 100 samples


        bpm = datafile.data('brachial_pressure_mean')
        i_zc = np.where((time > 60) & (bpm > bpm_base))[0]
        i_bpm = np.argmax(bpm[time > 60], axis = 0)
        
        
        if np.size(i_zc) == 0:
            # it did not fully compensate, lets take the maximum
            ws = s
            i_0 = i_bpm
            comp = False
        else:
            i_0 = i_zc[0]
            comp = True

        vol = datafile.data('addedVolume.volume')[i_0]/L2SI

        mean_rng = range(i_0 - int(np.floor(i_dt/2)), i_0 + int(np.floor(i_dt/2))) # use left and right helf interval for calculating a mean value
        gfr_m = np.mean(datafile.data('eGFR_m')[mean_rng])
        p_int_m = np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.p_isf')[mean_rng])/mmHg2SI
        vi_r =  np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.Vr')[mean_rng])
        vi_e = np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.V_excess')[mean_rng])/L2SI
        lymph_q = np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.lymph_flow')[mean_rng])/LpD2SI
        pcwp = np.mean(datafile.data('P_pv')[mean_rng])/mmHg2SI
        # edp = np.max()

        try:
            edv = np.max(datafile.data('EDV')[mean_rng])/L2SI*1000
            esv = np.min(datafile.data('ESV')[mean_rng])/L2SI*1000
        except KeyError:
            edv = 'NA'
            esv = 'NA'
        
        try:
            edp = np.max(datafile.data('EDP')[mean_rng])/mmHg2SI
            esp = np.min(datafile.data('ESP')[mean_rng])/mmHg2SI
        except KeyError:
            edp = 'NA'
            esp = 'NA'

        # aus funlib
        def calculateEF(volumes):
            esv = min(volumes)
            edv = max(volumes)
            return (edv - esv)/edv

        ef = calculateEF(datafile.data('V_LV')[i_0 - i_dt:i_0])*100
        co = datafile.data('CO')[i_0]*1000*60
        cvp = np.mean(datafile.data('P_sv')[i_0 - i_dt:i_0])/mmHg2SI
        hr = datafile.data('HR')[i_0]*60
        t = time[i_0]

        result = (i, vol, t, bpm[i_bpm], comp)
        result_set.append(result)

        ws = s.format(f_LV = i*imp_coeff, vol = vol,  time = t, eGFR = gfr_m, p_int = p_int_m, HR = hr, CO = co, CVP = cvp, PCWP = pcwp, EF = ef, EDV = edv, BPM = bpm[i_bpm]/mmHg2SI, Qlymph = lymph_q, Vint_r=vi_r, Vint_excess = vi_e,  ESP = esp, EDP = edp, comp = comp)

        file.write(ws)
        file.flush()
        print("Ok.")
        pass