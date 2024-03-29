# Gets the lowest cost from the OutputListingMain.txt and generates parametrized modelica model
import datetime
import matplotlib.pyplot as plt
import numpy

def readOutputListing(filename, runNumber = -1):
    """Reads the outputLIsting from GenOpt and returns params at lowest cost (or given run)
    
    Returns:
        (params, run, cost, datetime) 
    """

    with open(filename) as file:
        lines = file.readlines()

    param_names = {}
    run = -1
    cost = 1e300
    data_begin = False
    optim_date = 'generated at ' + str(datetime.datetime.now())

    var_set = {}
    var_set['cost'] = []
    for line in lines:
        
        if line.startswith('Start time:'):
            optim_date = 'optimized at ' + line[len('Start time: Thu '):]

        # prepare the header
        if not data_begin and line.startswith('Simulation Number'):
            # Simulation Number	Main Iteration	Step Number	f(x)	settings.phi0	settings.height	settings...
            cols = line.rstrip('\n').split('\t')
            for i in range(len(cols)):
                param_names[i] = cols[i]
                var_set[param_names[i]] = []
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
        var_set['cost'].append(cost)
        # list of ('name', value) tuples
        params = []
        
        # first 4 cols are runs and f(x) and last one is comment
        for col_num in range(4, len(cols)-1):
            param_name = param_names[col_num]
            # param_short_name = param_name[9:] if param_name.startswith('settings.') else param_name
            param_cur_val = float(cols[col_num])
            params.append((param_name, param_cur_val))
            var_set[param_name].append(param_cur_val)
        
        if run == runNumber:
            return (params, run, cost)
    
    return (params, run, cost, optim_date, var_set)
        

def generateModel(modelName, pack):
    """Generates modelice class parametrized by params
    
    Parameters:
        modelName   Base model
        pack        tuple (params, run, cost) where params is a list of (parameterName, value)
    """

    (params, run, cost, optim_date, var_set) = pack

    with open('%s_optimized.mo' % modelName, 'w') as file:
        header = """model {mn}_optimized "Generated by PostProcess/postprocess_optim.py {d} with cost {c} lowest at run {r}"
  extends {mn} (
    \n""".format(mn=modelName, d = optim_date, c = '%.6f' % cost, r = '%d' % run)
        footer = ');\nend %s_optimized;\n' % modelName

        param_lines = ('      %s = %.6e' % param for param in params)
        # use join to avoid placing the comma after the last
        lines = ',\n'.join(param_lines)


        file.write(header)
        file.writelines(lines)
        file.write(footer)        

def plotOptim(pack):
    (params, run, cost, optim_date, var_set) = pack
    fig, host = plt.subplots(figsize=(8,5)) # (width, height) in inches
    
    for (name, v) in var_set.items():
        if len(v) == 0:
            continue
        lb = numpy.min(v)
        ub = numpy.max(v)
        if lb == ub:
            v_r = list(0.5 for n in v)
        else:
            v_r = list((n-lb)/(ub - lb) for n in v)
        

        host.plot(v_r)

    plt.savefig('params_plot.png')


if __name__ is '__main__':
    pack = readOutputListing('OutputListingMain.txt')
    generateModel('SimpleExercise_ident', pack)
    plotOptim(pack)
    print("Done, JOne!")
    