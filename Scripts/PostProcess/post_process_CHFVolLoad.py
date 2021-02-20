import DyMat
import numpy as np
import time as timer
import os



filePattern = R"D:\Data\CVS_renalRegulation_CHF_baro_%0d.mat"
# filename = R"c:\home\UMICH\Valsalva\Results2\CardiovascularSystem.mat"


mmHg2SI = 133.322
ml2SI = 1e-6
bpm_base = 98*mmHg2SI
startTime = 60


s = '{f_LV}, {vol}, {time}, {eGFR}, {p_int}, {HR}, {CO}, {CVP}, {EF}, {BPM}, {comp}\n'
result_set = []


with open('CHF_VolumeLoading_baro.csv', 'w') as file:
    file.write(s.format(f_LV = 'f_LV', vol = 'vol', time = 'time', eGFR = 'eGFR', p_int = 'p_int', HR = 'HR', CO = 'CO', EF = 'EF', CVP = 'CVP', BPM = 'maxBP', comp = 'compensated'))

    for i in range(100, 0,-5):
        filename = filePattern % i
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
        i_dt = int(5/(time[100]/100)) # index of mean interval of 5s , measured over range of 100 samples


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

        vol = datafile.data('addedVolume.volume')[i_0]/ml2SI
        gfr_m = np.mean(datafile.data('eGFR_m')[i_0 - i_dt:i_0])
        p_int_m = np.mean(datafile.data('simplestLymphatic.p_int')[i_0 - i_dt:i_0])/mmHg2SI
        lymph_q = np.mean(datafile.data('simplestLymphatic.lymph_flow')[i_0 - i_dt:i_0])*1000*60*60*24

        # aus funlib
        def calculateEF(volumes):
            esv = min(volumes)
            edv = max(volumes)
            return (edv - esv)/edv

        ef = calculateEF(datafile.data('V_LV')[i_0 - i_dt:i_0])
        co = datafile.data('CO')[i_0]*1000*60
        cvp = np.mean(datafile.data('P_sv')[i_0 - i_dt:i_0])/mmHg2SI
        hr = datafile.data('HR')[i_0]*60
        t = time[i_0]

        result = (i, vol, t, bpm[i_bpm], comp)
        result_set.append(result)

        ws = s.format(f_LV = i, vol = vol,  time = t, eGFR = gfr_m, p_int = p_int_m, HR = hr, CO = co, CVP = cvp, EF = ef, BPM = bpm[i_bpm]/mmHg2SI, comp = comp)

        file.write(ws)
        print("Ok.")
        pass