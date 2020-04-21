import gen_dsin

dicSim = {}

params = 'Amref_factor, Vw_factor, ehm'
vars = params.split(',')
# strip the chars
vars = map(lambda x: x.strip(' ,()'), vars)
# prepare dictionary
# dicVars = {'Amref_factor': '%Amref_factor%'}
dicVars = {v : '%%%s%%' % v for v in vars}

gen_dsin.create_dsinORfinal(dicSim, dicVars, dsFileIn='dsin.txt', dsFileOut='dsinTemplate.txt')