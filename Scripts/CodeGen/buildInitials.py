""" Generates extending model with steady state parametrization, read from .mat simulation result.

Make sure all states listed in states.csv does not have initial equations.
Deleting state from the states.csv will disable generation of the initial value
Use writeInitStatesFromDsin() from ManipulateDsin to generate states.csv

MIT licence
Filip Ježek, University of Michigan. 2019
"""

import scipy.io as skipy
import DyMat
from matplotlib import pyplot as plt
import re
import TerminalDS
import ModelicaClass as mc
import os
import datetime

base_model_full_path = 'ADAN_main.SystemicTree.Variations.Renals.CVS_renalRegulation_HFrEF_15'
relative_folder = ''
exclude_filter = [] # ['Ra_phi', 'v_in', 'A']
# mid-cycle to avoid event collision during initialization 
steadyStateAt = 800
resetVPV = True
mat_file_path = 'settings.V_PV_init = 0.0021.mat' # or none to detect automatically

# get the main and path
if '.' in base_model_full_path:
    base_model_path, base_model = base_model_full_path.rsplit('.', 1)
else:
    base_model_path = None
    base_model = base_model_full_path


input_text = open('states.csv', 'r').read()
# incl_vars = 'volume|V_LV|V_RV|vol1|vol2|s|epsilon|fiSN|v_in'
incl_vars = '.+'
match = r'([\w.]+\.(?:' + incl_vars + r')),.+,(.+)'
m = re.findall(match, input_text)
lines = (line for line in m if line[0].rsplit('.', 1)[-1] not in exclude_filter)

# Build modelica Object Tree
mc_tree = mc.ModelicaClass.BuildObjectTree(lines, root=base_model)

if resetVPV:
    mc_tree.buildChildTree([['settings.V_PV_init', 'param']])
    vpvn = mc_tree.findNode('settings.V_PV_init')
    vpvn.start_val = 0

if mat_file_path is None:
    mat_file_path = relative_folder + base_model + '.mat'

mat_date = str(datetime.datetime.fromtimestamp(os.path.getmtime(mat_file_path)))

# open file with values
d = DyMat.DyMatFile(mat_file_path)
time = d.abscissa(2)[0]

nmsList = d.names(block = 2)
#nmsStr = "\r".join(nmsList)
steadyStateInd = TerminalDS.findLowestIndex(steadyStateAt, time, returnLast=True)

# pick value
for name in nmsList:
    # look for its existence in the states set
    n = mc_tree.findNode(name)
    if n is None:
        continue
    # if 'radial_vein_T3_R120' in name:
    #     print(d.data(name))
    n.start_val = d.data(name)[steadyStateInd]


mc_string = mc_tree.printObjectTree(indent_level=4)
init_model = base_model + '_init'
base_model_path_string = '' if base_model_path is None else base_model_path + '.'
print_string = \
    "model " + init_model + ' "Steady state initialization from ' + str(mat_date) + ' at time ' + str(time[steadyStateInd]) + '"\n' + \
    '  extends ' + base_model_path_string + mc_string + ';\n' + \
    'end ' + init_model + ';\n'

print(print_string)

    # write the output
with open(base_model + '_init.mo', 'w') as wf:
    wf.write(print_string)


print("Iam so DONE with this") 