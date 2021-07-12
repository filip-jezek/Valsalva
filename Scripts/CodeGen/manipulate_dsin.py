# -*- coding: utf-8 -*-
"""

@author: Scott Greenwood, Filip Jezek
This file has been grabbed from https://github.com/ORNL-Modelica/ModelicaPy/blob/master/modelicapy/gen_dsin.py at 460fa9d and modified by Filip Jezek, UMICH 2020
"""

import os
import re
from datetime import datetime
from typing import Iterable, Text

def readOutputListing(filename, runNumber = -1):
    """Reads the outputLIsting from GenOpt and returns params at lowest cost (or given run)

    Duplicated from from PostProcess/post_process_optim.py
    
    Returns:
        (params, run, cost, datetime) 
    """

    with open(filename) as file:
        lines = file.readlines()

    param_names = {}
    run = -1
    cost = 1e300
    data_begin = False
    optim_date = 'generated at ' + str(datetime.now())

    # line = []
    for line in lines:
        
        if line.startswith('Start time:'):
            optim_date = 'optimized at ' + line[len('Start time: Thu '):]

        # prepare the header
        if not data_begin and line.startswith('Simulation Number'):
            # Simulation Number	Main Iteration	Step Number	f(x)	settings.phi0	settings.height	settings...
            cols = line.rstrip('\n').split('\t')
            for i in range(len(cols)):
                param_names[i] = cols[i]
            data_begin = True
            continue
        elif not data_begin:
            # we are not yet there
            continue

        # find the parameters than are varied
        cols = line.rstrip('\n').split('\t')
        cur_cost = float(cols[3])

        if cur_cost > cost and run != runNumber:
            continue
        #else:

        # save this
        run = int(cols[0])
        
        cost = cur_cost
        # list of ('name', value) tuples
        params = []

        # first 4 cols are runs and f(x) and last one is comment
        for col_num in range(4, len(cols)-1):
            param_name = param_names[col_num]
            # param_short_name = param_name[9:] if param_name.startswith('settings.') else param_name
            param_cur_val = float(cols[col_num])
            params.append((param_name, param_cur_val))
        
        if run == runNumber:
            return (params, run, cost)
    
    return (params, run, cost, optim_date)

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
    """ Reads parameter initial value for all input parameters provided and returns a dict of (value, type)
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
                value_type_code = lin.split()[1]
                if value_type_code == '280':
                    value_type = 'stateInit'
                elif value_type_code == '272':
                    # a weird one, but still
                    value_type = 'stateInit2'
                elif value_type_code == '328':
                    value_type = 'discrete'
                elif value_type_code in ['361', '353']:
                    value_type = 'bool'
                elif value_type_code == '288':
                    value_type = 'guess'
                elif value_type_code == '264':
                    value_type = 'discreteInit'
                else:
                    value_type = value_type_code
                params[key] = (values[1], value_type)
                break
        if not isFound:
            raise ValueError('%s not found.' % key)

    return params

class OptimParam:
    """Just a simple data struct representing the optim parameter """
    
    def __init__(self, name, value = None, min = None, max = None, step = None, minmaxdev = 0.1):

        self.name = name
        self.value = value 
        self._min = min
        self._max = max
        self._step = step
        self.minmaxdev = minmaxdev
        # if value is not None:
        #     self.min = value*(1-minmaxdev) if min is None and self.value is not None else min
        #     self.max = value*(1+minmaxdev) if max is None else max
        #     self.step = value*0.01 if step is None else step
        # else:
        #     self.min = min
        #     self.max = max
        #     self.step = step

    @property
    def min(self):
        return self.value*(1-self.minmaxdev) if self._min is None and self.value is not None else self._min

    @property
    def max(self):
        return self.value*(1+self.minmaxdev) if self._max is None and self.value is not None else self._max

    @property
    def step(self):
        # return self._step
        # return self.value*0.01 if self._step is None and self.value is not None else self._step
        return 0 if self._step is None else self._step

def getInitParams(dsFileIn='dsin.txt', OptOutputFileIn = 'OutputListingMain.txt', paramsFile='dsin.txt') -> dict: 
    """ Gets dictionary of input parameters and its values

    returns dict of OptimParams
    """
    
    ops = {}
    if dsFileIn == paramsFile:
        with open(dsFileIn) as fil:
            lines_all = fil.readlines()

        lines_params = lines_all
        keys = getKeysfromDsin(lines_params)
        params = getValsFromDsin(lines_all, keys)    

        for k, p in (keys, params):
            ops[k] = OptimParam(k, p[0])

    else:
        with open(paramsFile) as fil:
            lines_params = fil.readlines()
            for l in lines_params:
                if l.startswith('#'):
                    continue
                s = l.rstrip('\n').split(',')
                if len(s) >= 5: 
                    op = OptimParam(s[0], min=s[2], max=s[3], step=s[4], value=s[1])                    
                elif len(s) == 4: 
                    op = OptimParam(s[0], min=s[2], max=s[3], value=s[1])
                elif len(s) == 2:
                    op = OptimParam(s[0], value=s[1])                    
                elif len(s) >= 1:
                    op = OptimParam(s[0])
                ops[s[0]] = op
        
        if dsFileIn is not None:
            with open(dsFileIn) as fil:
                lines_all = fil.readlines()
                params = getValsFromDsin(lines_all, ops.keys())
                for k, v in params.items():
                    ops[k].value = v[0]

        
        if OptOutputFileIn is not None:
            (params, _, _, _) = readOutputListing(OptOutputFileIn)
            for par in params:
                if par[0] not in ops:
                    print('WARINGN: ADDING PARAM %s!!' % par[0])
                    op = OptimParam(s[0])
                ops[par[0]].value = par[1]

    return ops

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

def getInputParams(dsFileIn='dsin.txt', filter = '', types = (280, 272), accept = (1, 2)):
    """ Search dsin for tunable Real parameters (type 280), optionally filtering by beggining

    Parameters

    filter  optional filter eg. filter='system.' returns params 'system.R' but not 'level.system.R'
    types   280 being Real parameter,  272 initial value
    accept  optional variable type, 1 being parameters and 2 for initial values only. Default is both.
    """

    tunable_params = []
    with open(dsFileIn) as file:
        lines = file.readlines()
        for line in lines:
            # seraches for pattern e.g. 'blbelbe 280   # settings.heart_R_LA[1]'
            res = re.search(r'(\d+)[ ]+(\d+)[ ]+# (%s[\w\.\[\]\d]+)' % filter, line)
            if res is not None:
                
                if int(res.groups()[1]) in types and int(res.groups()[0]) in accept:
                    tunable_params.append(res.groups()[2])
    
    return tunable_params

def writeInitStatesFromDsin(dsFileIn = 'dsin.txt', filter = '', outputFile = 'states.csv', accept = [2, 6], types = (280, 264, 272, 328, 361, 353)):
    initStates = getInputParams(dsFileIn, filter = filter, accept = accept, types = types)

    with open(dsFileIn) as file:
        lines = file.readlines()

    initStatesVals = getValsFromDsin(lines, initStates)
    initStatesLines = list('%s,%.3e,%s\n' % (k, v[0], v[1]) for (k, v) in initStatesVals.items())

    with open(outputFile, 'w') as file:
        file.write('# State variable;end value;\n')
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

def createDsinTemplate(keys, dsFileIn = 'dsin.txt', dsFileOut = 'dsinTemplate.txt', outputsOnly = False):
    template_mapping = getTemplateMapping(keys)
    
    # dicSim = {'StartTime': 0,'StopTime': 100}
    dicSim = {'levent': 0, 'evgrid' : 0}
    if outputsOnly:
        dicSim['lx'] = 0
        dicSim['lxd'] = 0
        dicSim['lxu'] = 0
        dicSim['lz'] = 0
        dicSim['lw'] = 0
        dicSim['la'] = 0
    # dicVars = {'m_flow': 10}

    create_dsinORfinal(dicSim, template_mapping, dsFileIn=dsFileIn, dsFileOut=dsFileOut)

def build_opt_command_file(filename, init_params:dict, step_frac = 0.05, run_type = 'sensitivity', ident_step_frac = 0.01):
    with open(filename, 'w') as file:
        file.write("""/* GenOpt command file 
Generated by manipulate_dsin.py at %s */

Vary{
""" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        for param_name, optimParam in init_params.items():
            value = float(optimParam.value)
            if run_type == 'sensitivity':
                step = '1'
            elif run_type == 'identification':
                step = optimParam.step if optimParam.step != 0 else value*(1 * ident_step_frac)

            min_ = optimParam.min if optimParam.min is not None else value*(1 - step_frac)
            max_ = optimParam.max if optimParam.max is not None else value*(1+ step_frac)
            
            line = '  Parameter{ Name = %s; Min = %s; Ini = %s; Max = %s; Step = %s; }\n' % (param_name, min_, value, max_ , step)
            file.write(line)

        if run_type == 'identification':
            file.write("""\n}
OptimizationSettings{
  MaxIte = 5000;
  MaxEqualResults = 4;
  WriteStepNumber = false;
  UnitsOfExecution = 4;
}\n""")
            if USEPSO:
                file.write("""\n
Algorithm{
  Main = GPSPSOCCHJ;
  NeighborhoodTopology = vonNeumann;
  NeighborhoodSize = 5;
  NumberOfParticle = 16;
  NumberOfGeneration = 16;
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

            else:
                file.write("""\n
 
Algorithm{
 Main = GPSHookeJeeves;
 MeshSizeDivider = 2;
 InitialMeshSizeExponent = 0;
 MeshSizeExponentIncrement = 1;
 NumberOfStepReduction = 6;
}
""" )
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

def writeInitParams(init_params:dict, paramsFile = 'init_params_default_vals.csv', step_frac_override = None, ident_frac_override = None):
    """ Writes params with default values into csv"""

    with open(paramsFile, 'w') as file:
        file.write('# name, value, [min, max, step]\n')
        for param, optimParam in init_params.items():

            if ident_frac_override is not None:
                min_ = optimParam.value - abs(optimParam.value)*ident_frac_override
            elif optimParam.min is None:
                # insert at least something
                min_ = optimParam.value - abs(optimParam.value)*0.5
            else:
                min_ = optimParam.min

            if ident_frac_override is not None:
                max_ = optimParam.value + abs(optimParam.value)*ident_frac_override
            elif optimParam.max is None:
                # insert at least something
                max_ = optimParam.value + abs(optimParam.value)*0.5
            else:
                max_ = optimParam.max

            if ident_frac_override is not None:
                step_ = optimParam.value*step_frac_override
            elif optimParam.step is None or optimParam.step == 0:
                # insert at least something
                step_ = optimParam.value*0.01
            else:
                step_ = optimParam.step

            s = '%s,%s,%s,%s,%s\n' % (param, optimParam.value, min_, max_, step_)
            file.write(s)

def prepareSA(paramsFile = 'params_for_SA.txt', regenerateParamsFromDsin = False, minMaxRange = 0):
    """

    params

    regenerateParamsFromDsin        rewrite the input file
    minMaxRange     fraction of minmax span. Keeps the loaded val if 0.    
    """
   
    # generate the params_for_SA.txt parameters list, which may be further edited.
    if regenerateParamsFromDsin:
        writeTunableParamsFromDsin(paramsFile)

    init_params = getInitParams(dsFileIn='dsin.txt', paramsFile=paramsFile)

    if minMaxRange > 0:
        for op in init_params.values():
            op._min = op.value*(1 + minMaxRange)
            op._max = op.value*(1 - minMaxRange)
    
    # writes the params with its initial value for simpler usage of other scripts, e.g. SA postprocessing
    writeInitParams(init_params, paramsFile = paramsFile)

    build_opt_command_file('opt_command_SA.txt', init_params, step_frac=0.05)

    createDsinTemplate(init_params, dsFileOut='dsinTemplate_SA.txt')

def prepareIdent(overrideFracs = False, regenerateParamsFromDsin = False, storeOnlyOutputs = False):
    paramsFile = 'params_for_ident.txt'
    # generate the params_for_SA.txt parameters list, which may be further edited.
    if regenerateParamsFromDsin:
        writeTunableParamsFromDsin(paramsFile)

    init_params = getInitParams(dsFileIn=DSFILEIN, OptOutputFileIn=OPTOUTPUTFILEIN, paramsFile=paramsFile)
    
    # writes the params with its initial value for simpler usage of other scripts, e.g. SA postprocessing
    if overrideFracs:
        writeInitParams(init_params, paramsFile=paramsFile, step_frac_override=0.1, ident_frac_override=0.5)
    else:
        writeInitParams(init_params, paramsFile=paramsFile)

    if overwriteOptParamFile:
        build_opt_command_file('opt_command.txt', init_params, run_type='identification')

    if overWriteDsinTemplate:
        createDsinTemplate(init_params, dsFileOut='dsinTemplate.txt', outputsOnly=storeOnlyOutputs)    
    

overwriteOptParamFile = True
overWriteDsinTemplate = True

# DSFILEIN = None
# OPTOUTPUTFILEIN = 'OutputListingMain.txt'

DSFILEIN = 'dsin.txt'
OPTOUTPUTFILEIN = None

USEPSO =  False

def run():
    # writeTunableParamsFromDsin('params_all.txt', filter='')
    # prepareSA(regenerateParamsFromDsin=False, minMaxRange=0.05)
    # prepareIdent(overrideFracs=False, regenerateParamsFromDsin=False, storeOnlyOutputs = False)
    writeInitStatesFromDsin(dsFileIn="dsin.txt")
    # writeTunableParamsFromDsin('params_all.txt', filter='')
    # writeTunableParamsFromDsin('params_settings.txt', filter='settings.')
    print('Done, Johne')
    
# if __name__ == "__main__":
run()
# writeInitStatesFromDsin(dsFileIn = 'dsin.txt', outputFile = 'params_for_SA.txt', filter = 'settings.', accept = [1, 2], types = (280, 272, 361))
    

