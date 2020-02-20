"""
Compares effects of tilt on the CVS

Creates dsinTilt00.txt - dsinTilt90.txt files and executes dymosim_Tilt_base model for different values of tilt and plots the results.
"""
import gen_dsin
import math
import matplotlib.pyplot as plt
import subprocess
import DyMat

def extractVars(filename):
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
    return var_set

regenerateInputFiles = False

if regenerateInputFiles:
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

# read the results
var_set00 = extractVars('dsres_Tilt00.mat')
var_set30 = extractVars('dsres_Tilt30.mat')
var_set60 = extractVars('dsres_Tilt60.mat')
var_set90 = extractVars('dsres_Tilt90.mat')

plot_var = 'Systemic1.brachial_L82_HeartLevel.u_C'
plot_labels = ('Supine', '30°', '60°', '90°')

fig1 = plt.figure()

plt.plot(var_set00['time'], var_set00[plot_var])
plt.plot(var_set30['time'], var_set30[plot_var])
plt.plot(var_set60['time'], var_set60[plot_var])
plt.plot(var_set90['time'], var_set90[plot_var])
plt.title(plot_var)
plt.legend(labels = plot_labels)

plt.show()