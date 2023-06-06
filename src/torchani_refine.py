#!/usr/bin/env python                                                                                                                                                                       # 2021-08-11

# Kiyoto Aramis Tanemura

# Once conformers are generated, refine them using ANI-2x potential using TorchANI
# Arg1: directory containing input mol file
# Arg2: directory for output xyz file

import os, sys
import json
from random import shuffle
#SDAS#
import  ase
from ase.io import read, write
from ase.optimize import BFGS, LBFGS 
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase import Atoms
#
from ase.calculators.calculator import Calculator, all_changes
#import ase.calculators.calculator
import ase.units
############
from ase.io import read
from ase.io.xyz import write_xyz
from ase.optimize import LBFGS
import torchani

if len(sys.argv) < 3:
    print('Refine conformer file by ANI-2x potential')
    print('Arg1: directory containing input mol file')
    print('Arg2: directory for output xyz file')
    quit()

inpath = sys.argv[1]
if inpath[-1] != '/':
    inpath += '/'

outpath = sys.argv[2]
if outpath[-1] != '/':
    outpath += '/'

filelist = [x for x in os.listdir(inpath) if x[-3:] == 'mol']
shuffle(filelist)

#calculator = torchani.models.ANI2x().ase()
##SDAS
calculator = torchani.models.ANI2x(periodic_table_index=True).double()
#geometry.set_calculator(model.ase())
###
for thefile in filelist:
    outfile = thefile.split('.')[0] + '.xyz'
    if outfile in os.listdir(outpath):
        continue
    ase_atoms = read(inpath + thefile)
   # ase_atoms.set_calculator(calculator)
    ase_atoms.set_calculator(calculator.ase())
    dyn = LBFGS(ase_atoms)
    dyn.run(fmax=0.001, steps = 2000)
    with open(outpath + outfile, 'w') as f:
        write_xyz(f, [ase_atoms])
#    ase_atoms.write(outpath + outfile, format = 'xyz', plain=True)
