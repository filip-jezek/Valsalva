import DyMat
import numpy as np


filePattern = R"D:\Data\CVS_renalRegulation_CHF_%0d.mat"
filename = R"c:\home\UMICH\Valsalva\Results2\CardiovascularSystem.mat"


mmHg2SI = 133.322
bpm_base = 92*mmHg2SI
startTime = 60


s = '{f_LV}, {vol}, {time}, {BPM}, {comp}'
result_set = []

with open('VolumeLoadingLog', 'w') as file:
    file.write(s.format(f_LV = 'f_LV', vol = 'vol', time = 'time', BPM = 'maxBP', comp = 'compensated'))

    for i in range(100, 0,-10):
        filename = filePattern % i
        filename = R"c:\home\UMICH\Valsalva\Results2\CardiovascularSystem.mat"
    
        try:
            datafile = DyMat.DyMatFile(filename)
        except FileNotFoundError:
            print("File %0d  not found" % i)
            continue


        time = datafile.abscissa(2)[0]

        bpm = datafile.data('brachial_pressure_mean')
        i_zc = np.where((time > 60) & (bpm > bpm_base))
        i_bpm = np.argmax(bpm[time > 60], axis = 0)
        
        
        if np.size(i) == 0:
            ws = s
            i_0 = i_bpm
            comp = False
        else:
            i_0 = i_zc[0][0]
            comp = True

        vol = datafile.data('addedVolume.volume')[i_0]
        # vol = datafile.data('totalVolume')[i_0]*1000
        t = time[i_0]

        result = (i, vol, t, bpm[i_bpm], comp)
        result_set.append(result)

        ws = s.format(f_LV = i, vol = vol,  time = t, BPM = bpm[i_bpm])

        file.write(ws)
        
        pass




