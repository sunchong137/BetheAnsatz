import sys
import numpy as np
sys.path.append('../')
from BAFT import *
import os

U = float(sys.argv[1])
write2file = False
outdir = "data/"
if not os.path.exists(outdir):
    os.mkdir(ourdir)
beta = np.array([0.01, 0.001, 0.0001])
#beta = np.linspace(0.1,10,100,endpoint=True)
Tgrid = 1./beta
#Tgrid = np.linspace(0.00,2.0,41,endpoint=True)
#Tgrid = np.linspace(1.05,2.0,20, endpoint=True)
solve_energy_curve(U,Tgrid,savefile=write2file)


