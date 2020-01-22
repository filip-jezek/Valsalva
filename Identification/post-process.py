import scipy.io as skipy
import DyMat
from matplotlib import pyplot as plt
import re
import TerminalDS



# filename = 'CardiovascularSystem_ADANExport'

steadyStateAt = 30
filename = 'dsres.mat'
# mat = skipy.loadmat(filename)
d = DyMat.DyMatFile(filename)
time = d.abscissa(2)[0]

with open('dsout.out', 'w') as file:
	file.write("banik pico")
