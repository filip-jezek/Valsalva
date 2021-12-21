import DyMat
import numpy as np
import time as timer
import os


mmHg2SI = 133.322
L2SI = 1e-3
LpD2SI = 1/(1000*60*60*24)
vol_norm = 3.697 # L

hftypes = ['CVS_Ols_nobaro_LV100_']
# hftypes = ['CVS_nobaro_Lumens_LV100_', 'CVS_nobaro_Lumens_LV50_', 'CVS_nobaro_Lumens_LV25_']

for hftype in hftypes:
    filenames = [
        'vol4.mat',
        'vol3.6.mat',
        'vol3.2.mat',
        'vol2.8.mat',
        'vol2.4.mat',
        'vol2.mat',
        'vol1.6.mat',
        'vol1.2.mat',
        'vol0.8.mat',
        'vol0.4.mat',
        'vol0.mat',
        'vol-0.4.mat',
        'vol-0.8.mat',
        'vol-1.2.mat',
        'vol-1.6.mat',
        'vol-2.mat',]

    exp_range = range(0, len(filenames))

    s = '{f_LV}, {vol}, {time}, {eGFR}, {p_int}, {HR}, {CO}, {CVP}, {PCWP}, {EF}, {EDV}, {BPm}, {BPs}, {BPd}, {Qlymph}, {Vint_r},{Vint_excess}, {ESP}, {EDP}, {Pwr_RV}, {Pwr_LV}, {dCar_d}, {dAor_d}\n'
    result_set = []

    outputFile = 'VolLoad_result_%s.csv' % hftype
    with open(outputFile, 'w') as file:
    # with open('CHF_VolumeLoading.csv', 'w') as file:    
        file.write(s.format(f_LV = 'f_LV', vol = 'vol', time = 'time', eGFR = 'eGFR', p_int = 'p_int', HR = 'HR', CO = 'CO', EF = 'EF', EDV = 'EDV', CVP = 'CVP', PCWP = 'PCWP', BPm = 'BPm', BPs = 'BPs', BPd = 'BPd', Qlymph = 'Qlymph', Vint_r= 'V_Isf_Rel', Vint_excess = 'V_Isf_Exc', ESP = 'ESP', EDP = 'EDP', Pwr_LV = 'Pwr_LV', Pwr_RV = 'Pwr_RV', dCar_d = 'dCar_d', dAor_d = 'dAor_d'))

        for i in exp_range:
            # if hftype == '':
            filename = hftype + filenames[i]
            # else:
            #     filename = '\\'.join([hftype, filenames[i]])
            # filename = R"c:\home\UMICH\Valsalva\Results2\CardiovascularSystem.mat"

            try:
                print("Loading file %s..." % filename, end = '')
                tic = timer.time()
                datafile = DyMat.DyMatFile(filename)
                if datafile.abscissa(2)[0][-1] < 10:
                    raise IndexError()
                print("Ok [%0.3fs]. Processing..." % (timer.time() - tic), end = '')                
            except FileNotFoundError:
                print(" File %0d  not found !" % i)
                continue
            except:
                print(" File %0d  failed, show must go on" % i)

                # os.system("c:\\Program Files\\Dymola 2021\\bin\\dsres2dsf %s %s" % (filename, filename.replace('.mat', '.sdf'))

                continue


            time = datafile.abscissa(2)[0]
            i_dt = int(6/(time[100]/100)) # index of mean interval of 6s , measured over range of 100 samples

            # from funlib
            def findLowestIndex(time, timeArr):
                lst = timeArr.tolist()
                return next((i for i, x in enumerate(lst) if x >= time))

            # from funlib
            def findInterval(t_from, t_to, timeArr):
                return range(findLowestIndex(t_from, timeArr), findLowestIndex(t_to, timeArr))

            mean_rng = findInterval(time[-1] - 4, time[-1] -1, time)

            # imp_coeff = datafile.data('settings.V_PV_init')[0]
            
            vol = np.mean(datafile.data('totalVolume')[mean_rng])/L2SI
            imp_coeff = vol - vol_norm

            gfr_m = np.mean(datafile.data('eGFR_m')[mean_rng])
            p_int_m = np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.p_isf')[mean_rng])/mmHg2SI
            vi_r =  np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.Vr')[mean_rng])
            vi_e = np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.V_excess')[mean_rng])/L2SI
            lymph_q = np.mean(datafile.data('simplestLymphaticDynamicSpeedUp.lymph_flow')[mean_rng])/LpD2SI
            pcwp = np.mean(datafile.data('P_pv')[mean_rng])/mmHg2SI
            BPm = np.mean(datafile.data('brachial_pressure_mean')[mean_rng])/mmHg2SI
            BPs = np.max(datafile.data('brachial_pressure')[mean_rng])/mmHg2SI
            BPd = np.min(datafile.data('brachial_pressure')[mean_rng])/mmHg2SI
            Pwr_LV = np.max(datafile.data('heartComponent.ventricles.power_LV')[mean_rng])
            Pwr_RV = np.max(datafile.data('heartComponent.ventricles.power_RV')[mean_rng])
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

            ef = calculateEF(datafile.data('V_LV')[mean_rng])*100
            co = np.mean(datafile.data('CO')[mean_rng])*1000*60
            cvp = np.mean(datafile.data('P_sv')[mean_rng])/mmHg2SI
            hr = np.mean(datafile.data('HR')[mean_rng])*60
            dcar_d = np.max(datafile.data('SystemicComponent.baroreflex_system.carotid_distention')[mean_rng]) - np.min(datafile.data('SystemicComponent.baroreflex_system.carotid_distention')[mean_rng])
            daor_d = np.max(datafile.data('SystemicComponent.baroreflex_system.aortic_distention')[mean_rng]) - np.min(datafile.data('SystemicComponent.baroreflex_system.aortic_distention')[mean_rng])
            
            t = time[-1]

            ws = s.format(f_LV = imp_coeff, vol = vol,  time = t, eGFR = gfr_m, p_int = p_int_m, HR = hr, CO = co, CVP = cvp, PCWP = pcwp, EF = ef, EDV = edv, BPm = BPm, BPs = BPs, BPd = BPd, Qlymph = lymph_q, Vint_r=vi_r, Vint_excess = vi_e,  ESP = esp, EDP = edp, Pwr_LV = Pwr_LV, Pwr_RV = Pwr_RV, dCar_d = dcar_d, dAor_d = daor_d)

            file.write(ws)
            file.flush()
            print("Ok.")
    pass
print('Thats it.')