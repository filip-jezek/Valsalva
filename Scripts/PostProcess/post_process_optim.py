# Gets the lowest cost from the OutputListingMain.txt and generates parametrized modelica model
import datetime

def readOutputListing(filename, runNumber = -1):
    """Reads the outputLIsting from GenOpt and returns params at lowest cost (or given run)
    
    Returns:
        (params, run, cost) 
    """

    with open(filename) as file:
        lines = file.readlines()

    param_names = {}
    run = -1
    cost = 1e300
    data_begin = False
    # line = []
    for line in lines:

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
    
    return (params, run, cost)
        

def generateModel(modelName, pack):
    """Generates modelice class parametrized by params
    
    Parameters:
        modelName   Base model
        pack        tuple (params, run, cost) where params is a list of (parameterName, value)
    """

    (params, run, cost) = pack

    with open('%s_optimized.mo' % modelName, 'w') as file:
        header = """model {mn}_optimized "Generated by PostProcess/postprocess_optim.py at {d} with cost {c} lowest at run {r}"
  extends {mn} (\n""".format(mn=modelName, d = datetime.datetime.now(), c = '%.6f' % cost, r = '%d' % run)
        footer = ');\nend %s;\n' % modelName

        param_lines = ('    %s = %.6e' % param for param in params)
        # use join to avoid placing the comma after the last
        lines = ',\n'.join(param_lines)


        file.write(header)
        file.writelines(lines)
        file.write(footer)


if __name__ is '__main__':
    pack = readOutputListing('OutputListingMain_optim.txt')
    generateModel('OlufsenTriSeg', pack)