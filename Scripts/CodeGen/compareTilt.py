"""
Compares effects of tilt on the CVS and plots the results

Creates dsinTilt00.txt - dsinTilt90.txt files and executes dymosim_Tilt_base model for different values of tilt and plots the results.
"""
import gen_dsin
import math
import matplotlib.pyplot as plt
import subprocess
import DyMat
import time

def extractVars(filename, tic = time.time()):
    # use just a subset or use all
    # vrs_names = ['Systemic1.posterior_tibial_T4_R236.u_C', 'Systemic1.aortic_arch_C2.port_a.pressure']
    # vrs_names = d.names(block = 2)
    d = DyMat.DyMatFile(filename)
    vrs_names = d.names()
    timevar = d.abscissa(2)[0]

    var_set = {}
    var_set['time'] = timevar
    for vn in vrs_names:
        var_set[vn] = d.data(vn)
    
    print("Opening result ", filename, " in ", time.time() - tic, " s")
    return var_set

def generateInputFiles():
    # create Tilt 0
    gen_dsin.create_dsinORfinal({}, {'Tilt_ramp.height': 0}, dsFileIn='dsinTilt.txt', dsFileOut='dsinTilt00.txt')
    # create Tilt 30
    gen_dsin.create_dsinORfinal({}, {'Tilt_ramp.height': math.pi/6}, dsFileIn='dsinTilt.txt', dsFileOut='dsinTilt30.txt')
    # create Tilt 60
    gen_dsin.create_dsinORfinal({}, {'Tilt_ramp.height': math.pi/3}, dsFileIn='dsinTilt.txt', dsFileOut='dsinTilt60.txt')
    # create Tilt 90
    gen_dsin.create_dsinORfinal({}, {'Tilt_ramp.height': math.pi/2}, dsFileIn='dsinTilt.txt', dsFileOut='dsinTilt90.txt')

    # run dymosims
    subprocess.run(['dymosim_Tilt_base.exe', 'dsinTilt00.txt', 'dsres_Tilt00.mat'], shell=True, check = True)
    subprocess.run(['dymosim_Tilt_base.exe', 'dsinTilt30.txt', 'dsres_Tilt30.mat'], shell=True, check = True)
    subprocess.run(['dymosim_Tilt_base.exe', 'dsinTilt60.txt', 'dsres_Tilt60.mat'], shell=True, check = True)
    subprocess.run(['dymosim_Tilt_base.exe', 'dsinTilt90.txt', 'dsres_Tilt90.mat'], shell=True, check = True)

def plot_tiltvars(var_sets, plot_var, plot_labels = None, plot_title = None, ylim = None, subplot = None):
    if plot_title is None:
        plot_title = plot_var

    (var_set00, var_set30, var_set60, var_set90) = var_sets

    if subplot is None:
        plt.figure()
    else:
         plt.subplot(subplot)

    plt.plot(var_set00['time'], var_set00[plot_var]/133.32)
    plt.plot(var_set30['time'], var_set30[plot_var]/133.32)
    plt.plot(var_set60['time'], var_set60[plot_var]/133.32)
    plt.plot(var_set90['time'], var_set90[plot_var]/133.32)
    plt.title(plot_var)
    if plot_labels is not None:
        plt.legend(labels = plot_labels, ncol = len(plot_labels))
    plt.title(plot_title)
    plt.xlim(28.3, 29.5)
    if ylim is not None:
        plt.ylim(ylim)
    

def plot_tiltpos(var_set, plot_vars, legend = None, plot_title = None, ylim = None, subplot = None):
    
    if subplot is None:
        plt.figure()
    else:
         plt.subplot(subplot)

    for plot_var in plot_vars:
        plt.plot(var_set['time'], var_set[plot_var]/133.32)
    plt.title(plot_var)
    if legend is not None:
        plt.legend(labels = legend, loc = 'right')
    plt.title(plot_title)
    plt.xlim(28.3, 29.5)
    if ylim is not None:
        plt.ylim(ylim)
    

regenerateInputFiles = False
tic = time.time()

if regenerateInputFiles:
    generateInputFiles()


# read the results
var_set00 = extractVars('dsres_Tilt00.mat')
var_set30 = extractVars('dsres_Tilt30.mat')
var_set60 = extractVars('dsres_Tilt60.mat')
var_set90 = extractVars('dsres_Tilt90.mat')

plt.close('All')

fig1 = plt.figure(1)
plt.clf()
var_sets = (var_set00, var_set30, var_set60, var_set90)
plot_labels = ('0°', '30°', '60°', '90°')
plot_tiltvars(var_sets, 'Systemic1.brachial_L82_HeartLevel.u_C', plot_labels = plot_labels, plot_title = 'brachial_L82_HeartLevel', ylim=(60, 140), subplot = 321)
plot_tiltvars(var_sets, 'Systemic1.femoral_R226.u_C', plot_labels = plot_labels, plot_title = 'femoral_R226', ylim=(60, 180), subplot = 322)
plot_tiltvars(var_sets, 'Systemic1.internal_carotid_R8_B.u_C', plot_labels = plot_labels, plot_title = 'internal_carotid_R8_B', ylim=(60, 200), subplot = 323)
plot_tiltvars(var_sets, 'Systemic1.thoracic_aorta_C108.u_C', plot_labels = plot_labels, plot_title = 'thoracic_aorta_C108', ylim=(60, 200), subplot = 324)
plot_tiltvars(var_sets, 'Systemic1.tibiofibular_trunk_R234.u_C', plot_labels = plot_labels, plot_title = 'tibiofibular_trunk_R234', ylim=(60, 200), subplot = 325)
plot_tiltvars(var_sets, 'heartComponent.aorticValve.q_out.pressure', plot_labels = plot_labels, plot_title = 'aorticValve', ylim=(60, 200), subplot = 326)
plt.tight_layout()

fig2 = plt.figure(2)
plt.clf()
plot_vars = ('Systemic1.brachial_L82_HeartLevel.u_C',
'Systemic1.femoral_R226.u_C',
'Systemic1.internal_carotid_R8_B.u_C',
'Systemic1.thoracic_aorta_C108.u_C',
'Systemic1.tibiofibular_trunk_R234.u_C',
'heartComponent.aorticValve.q_out.pressure')
legend = ('brachial_L82_HL',
'femoral_R226',
'internal_carotid_R8_B',
'thoracic_aorta_C108',
'tibiofibular_trunk_R234',
'aorticValve')
plot_titles = ('Supine', '30°', '60°', '90°')

plot_tiltpos(var_set00, plot_vars, plot_title = plot_titles[0], subplot = 221)
plot_tiltpos(var_set30, plot_vars, legend = legend, plot_title = plot_titles[1], subplot = 222)
plot_tiltpos(var_set60, plot_vars, plot_title = plot_titles[2], subplot = 223)
plot_tiltpos(var_set90, plot_vars, plot_title = plot_titles[3], subplot = 224)
plt.tight_layout()

plot_vars_position = ('Systemic1.brachial_L82_HeartLevel.hydrostatic_level',
'Systemic1.femoral_R226.hydrostatic_level',
'Systemic1.internal_carotid_R8_B.hydrostatic_level',
'Systemic1.thoracic_aorta_C108.hydrostatic_level',
'Systemic1.tibiofibular_trunk_R234.hydrostatic_level',
'None')

for i in range(len(var_sets)):
    s = print('For %s tilt:' % plot_labels[i])
    for j in range(len(plot_vars)-1):
        print("  %s at %s cm" % (legend[j], str(round(var_sets[i][plot_vars_position[j]][-1]*100, 1)) ))
        
plt.show()