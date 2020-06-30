# -*- coding: utf-8 -*-
"""

@author: Scott Greenwood, Filip Jezek
This file has been grabbed from https://github.com/ORNL-Modelica/ModelicaPy/blob/master/modelicapy/gen_dsin.py at 460fa9d and modified by Filip Jezek, UMICH 2020
"""

import os
import re
from datetime import datetime
from typing import Iterable, Text

def getKeysfromDsin(lines):
    """ Reads initial parameter names from dymola dsin file or from any text file providing one param per line
    """

    if lines[0].startswith('#1'):
        # marker for dsin type of file, so we would have to skip initial gibberish
        initNamesSegment = False
    else:
        # lets start right over
        initNamesSegment = True

    keys = []
    for i in range(len(lines)):
        if not initNamesSegment and lines[i].startswith('char initialName'):
            # marks begining of the segment
            initNamesSegment = True
            continue
        if initNamesSegment and (lines[i] == '\n' or lines[i].startswith('double initialValue')):
            # end of the segment
            break
        if initNamesSegment:
            # add to list but trim the line breaks
            keys.append(lines[i].rstrip('\n'))
    
    return keys

def getValsFromDsin(lines:Iterable[Text], keys: Iterable[Text]):
    """ Reads parameter initial value for all input parameters provided
    """

    params = {}
    for key in keys:
        isFound = False
        for x in range(len(lines)):
            lin = lines[x]
            # we have to mark the key so its not found anywhere else
            marked_key = '# ' + key
            if marked_key + '\n' in lin:
                isFound = True
                value_line = lines[x-1]
                values = list(map(float, value_line.split()))
                params[key] = values[1]
                break
        if not isFound:
            raise ValueError('%s not found.' % key)

    return params

def getInitParams(dsFileIn='dsin.txt', paramsFile='dsin.txt') -> dict: 
    """ Gets dictionary of input parameters and its values
    """
    with open(dsFileIn) as fil:
        lines_all = fil.readlines()

    if dsFileIn == paramsFile:
        lines_params = lines_all
    else:
        with open(paramsFile) as fil:
            lines_params = fil.readlines()

    keys = getKeysfromDsin(lines_params)
    params = getValsFromDsin(lines_all, keys)
    return params

def create_dsinORfinal(dicSim, dicVars, dsFileIn='dsin.txt', dsFileOut='dsinTemp.txt', simulator='dymola', autoRewrite = True):
    '''
    Generate a new dsin.txt file with modified values from discSim and dicVars
    from a dsin.txt or dsfinal.txt file with a name specified by dsFileOut.

    dsFileIn = dsfile to be read (e.g., dsin.txt)
    simulator = currently only supports dymola
    dsFileOut = dsfile to be written
    dicSim = simulation settings to be changed
    dicVars = inital variable values to be changed

    Example:
        create_dsinORfinal({'StartTime': 0,'StopTime': 100},
                     {'var1.k': 10},
                     'dsfinal.txt',
                     'dsin.txt',
                     'dymola')                    
    '''

    if simulator != "dymola":
        raise ValueError('Argument "simulator" needs to be set to "dymola".')

    if not os.path.isfile(dsFileIn):
        raise IOError("File {} does not exist".format(dsFileIn))

    if dsFileIn == dsFileOut and autoRewrite == False:
        answer = input('''Input and ouput file names match.
    The input file will be overwritten. Continue [y/n]? ''')
        if answer.lower() in ['y', 'yes']:
            # Do nothing
            pass
        elif answer.lower() in ['n', 'no']:
            print('Program terminated')
            return
        else:
            raise IOError('Response not recognized. Program terminated')

    # Add '#' identifier used in dsin/dsfinal files to dicVars_i
    
    dicVars_i = {}
    for key, value in dicVars.items():
        dicVars_i['# %s' % key] = value

    with open(dsFileIn) as fil:
        lines = fil.readlines()

    # Dymola generated dsin.txt files have an extra blank line
    # near the top that can be used to differentiate it from dsfinal.txt
    if len(lines[6]) == 1:
        # dsin.txt file type
        dstype = 0
    else:
        # dsfinal.txt file type
        dstype = 1

    # Experiment parameter options:
    opts_Exp = {'StartTime': (9, 8), 'StopTime': (11, 9),
                'Increment': (12, 10), 'nInterval': (13, 11),
                'Tolerance': (14, 12), 'MaxFixedStep': (16, 13),
                'Algorithm': (18, 14)}

    # Method tuning parameter options:
    opts_Tun = {'grid': (44, 18), 'nt': (54, 19),
                'dense': (55, 20), 'evgrid': (56, 21),
                'evu': (57, 22), 'evuord': (58, 23),
                'error': (59, 24), 'jac': (60, 25),
                'xd0c': (61, 26), 'f3': (62, 27),
                'f4': (63, 28), 'f5': (64, 29),
                'debug': (65, 30), 'pdebug': (66, 31),
                'fmax': (67, 32), 'ordmax': (68, 33),
                'hmax': (69, 34), 'hmin': (70, 35),
                'h0': (71, 36), 'teps': (72, 37),
                'eveps': (73, 38), 'eviter': (74, 39),
                'delaym': (75, 40), 'fexcep': (76, 41),
                'tscale': (77, 42)}

    #  Output parameter options:
    opts_Out = {'lprec': (87, 48), 'lx': (88, 49),
                'lxd': (89, 50), 'lu': (90, 51),
                'ly': (91, 52), 'lz': (92, 53),
                'lw': (93, 54), 'la': (94, 55),
                'lperf': (95, 56), 'levent': (96, 57),
                'lres': (97, 58), 'lshare': (98, 59),
                'lform': (99, 60)}

    fnew = open(dsFileOut, 'w')
    # Simulation Settings:
    # Search for each of the experiment options and replace the value
    for key, value in opts_Exp.items():
        if key in dicSim:
            lines[value[dstype]] = ' %s\n' % dicSim[key]

    # Search for each of the method tuning options and replace the value
    for key, value in opts_Tun.items():
        if key in dicSim:
            lines[value[dstype]] = ' %s\n' % dicSim[key]

    # Search for each of the output options and replace the value
    for key, value in opts_Out.items():
        if key in dicSim:
            lines[value[dstype]] = ' %s\n' % dicSim[key]
    # End Simulation Settings:

    # Initial Condition Settings:
    # Search for variable and replace the inital value
    # Variables are defined by 6 columns. Only 4 are currently read
    # Guidance is taken from section starting with:
    # "Matrix with 6 columns defining the initial value calculation"
    # c1 : type of initial value
    # c2 : value
    # c3 : Minimum value (ignored, if Minimum >= Maximum)
    # c4 : Maximum value (ignored, if Minimum >= Maximum)
    # c5 : Category of variable
    # c6 : Data type of variable and flags according to dsBaseType
    for key, value in dicVars_i.items():
        isFound = False
        for x in range(len(lines)):
            lin = lines[x]
            # we have to look up to end of line, due to parameters with the same prefix, e.g. amp_factor and amp_factor_min
            if key + '\n' in lin:
                isFound = True
                temp = lines[x-1]
                temp = list(map(float, temp.split()))
                c1 = int(temp[0])
                c2 = dicVars_i[key]
                c3 = temp[2]
                c4 = temp[3]
                lines[x-1] = ' %s %s %s %s\n' % (c1, c2, c3, c4)
        if not isFound:
            print(key)
            print(lin)
            raise ValueError('%s not found in %s.' % (key, dsFileIn))

    fnew.writelines(lines)
    fnew.close()

def getInputParams(dsFileIn='dsin.txt', filter = '', accept = (1, 2)):
    """ Search dsin for tunable Real parameters (type 280), optionally filtering by beggining

    Parameters

    filter  optional filter eg. filter='system.' returns params 'system.R' but not 'level.system.R'
    accept  optional variable type, 1 being parameters and 2 for initial values only. Default is both.
    """

    tunable_params = []
    with open(dsFileIn) as file:
        lines = file.readlines()
        for line in lines:
            # seraches for pattern e.g. 'blbelbe 280   # settings.heart_R_LA'
            res = re.search(r'(\d+)[ ]+(\d+)[ ]+# (%s[\w\.]+)' % filter, line)
            if res is not None:
                
                if res.groups()[1] == '280' and int(res.groups()[0]) in accept:
                    tunable_params.append(res.groups()[2])
    
    return tunable_params

def writeInitStatesFromDsin(dsFileIn = 'dsin.txt', outputFile = 'states.csv'):
    initStates = getInputParams(dsFileIn, accept = [2])

    with open(dsFileIn) as file:
        lines = file.readlines()

    initStatesVals = getValsFromDsin(lines, initStates)
    initStatesLines = list('%s;%.3e;\n' % (k, v) for (k, v) in initStatesVals.items())

    with open(outputFile, 'w') as file:
        file.write('State variable;end value;\n')
        file.writelines(initStatesLines)

def writeTunableParamsFromDsin(outputFile, filter='settings.'):
    " Out of dsin.txt generates list of tunable parameters. Usually enough to run once."
    tunable_params = getInputParams(filter = filter, accept = [1])
    with open(outputFile, 'w') as file:
        file.writelines('\n'.join(tunable_params))


# prepare the dsin for genopt
def getTemplateMapping(keys: dict) ->dict:
    """Gets the mapping key: '%key%'
    """
    return {key:'%%%s%%' % key for key in keys.keys()}

def createDsinTemplate(keys, dsFileIn = 'dsin.txt', dsFileOut = 'dsinTemplate.txt'):
    template_mapping = getTemplateMapping(keys)
    
    # dicSim = {'StartTime': 0,'StopTime': 100}
    dicSim = {}
    # dicVars = {'m_flow': 10}

    create_dsinORfinal(dicSim, template_mapping, dsFileIn=dsFileIn, dsFileOut=dsFileOut)

def build_opt_command_file(filename, init_params:dict, step_frac = 0.05, run_type = 'sensitivity', ident_step_frac = 0.01):
    with open(filename, 'w') as file:
        file.write("""/* GenOpt command file 
Generated by manipulate_dsin.py at %s */

Vary{
""" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        for param_name, param_value in init_params.items():
            value = float(param_value)
            if run_type == 'sensitivity':
                inputs = (param_name, value*(1 - step_frac), value, value*(1+ step_frac), '1')
            elif run_type == 'identification':
                inputs = (param_name, value*(1 - step_frac), value, value*(1+ step_frac), '%.3e' % (value*ident_step_frac))

            line = '  Parameter{ Name = %s; Min = %s; Ini = %s; Max = %s; Step = %s; }\n' % inputs
            file.write(line)

        if run_type == 'identification':
            file.write("""\n}
OptimizationSettings{
  MaxIte = 5000;
  MaxEqualResults = 4;
  WriteStepNumber = false;
  UnitsOfExecution = 4;
}

Algorithm{
  Main = GPSPSOCCHJ;
  NeighborhoodTopology = vonNeumann;
  NeighborhoodSize = 5;
  NumberOfParticle = 16;
  NumberOfGeneration = 20;
  Seed = 1;
  CognitiveAcceleration = 2.8;
  SocialAcceleration = 1.3;
  MaxVelocityGainContinuous = 0.5;
  MaxVelocityDiscrete = 4;
  ConstrictionGain = 0.5;
  MeshSizeDivider = 2;
  InitialMeshSizeExponent = 0;
  MeshSizeExponentIncrement = 1;
  NumberOfStepReduction = 4;
}""" )
            # jsut to get back at proper  indent
            pass

        else:        
            file.write("""\n}

OptimizationSettings{
MaxIte = 500;
MaxEqualResults = 10;
WriteStepNumber = false;
UnitsOfExecution = 4;
}

Algorithm {
Main = Parametric;
StopAtError = false;
}
""")

    # jsut to get back at proper indent at the method end
    pass

def writeInitParams(init_params:dict):
    """ Writes params with default values into csv"""

    with open('init_params_default_vals.csv', 'w') as file:
        for param, val in init_params.items():
            s = '%s,%s\n' % (param, val)
            file.write(s)

def prepareInitStates():
    writeInitStatesFromDsin()

def prepareSA():
    paramsFile = 'params_for_SA.txt'
    # generate the params_for_SA.txt parameters list, which may be further edited.
    # uncomment if thats the first run to get all tunable parameters
    # writeTunableParamsFromDsin(paramsFile)

    init_params = getInitParams(dsFileIn='dsin.txt', paramsFile=paramsFile)
    
    # writes the params with its initial value for simpler usage of other scripts, e.g. SA postprocessing
    writeInitParams(init_params)

    build_opt_command_file('opt_command_SA.txt', init_params, step_frac=0.05)

    createDsinTemplate(init_params, dsFileOut='dsinTemplate_SA.txt')

def prepareIdent():
    paramsFile = 'params_for_ident.txt'
    # generate the params_for_SA.txt parameters list, which may be further edited.
    # uncomment if thats the first run to get all tunable parameters
    # writeTunableParamsFromDsin(paramsFile)

    init_params = getInitParams(dsFileIn='dsin.txt', paramsFile=paramsFile)
    
    # writes the params with its initial value for simpler usage of other scripts, e.g. SA postprocessing
    writeInitParams(init_params)

    build_opt_command_file('opt_command.txt', init_params, step_frac=0.5, run_type='identification', ident_step_frac=0.01)

    createDsinTemplate(init_params, dsFileOut='dsinTemplate.txt')    

if __name__ == "__main__":

    prepareSA()
    # prepareIdent()
    print('Done, Johne')

