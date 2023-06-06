#!/usr/bin/env python                                                                                                                                                                                

# 2021-04-06

# Kiyoto Aramis Tanemura

# Once conformers are generated, refine them using ANI-2x potential.
# Arg1: directory containing input mol file
# Arg2: directory for output xyz file

import os, sys
import json
from random import shuffle
from ase.io import read
from ase.optimize import LBFGS
with open('paths.json', 'r') as f:
    paths = json.load(f)
if not os.path.isdir(paths['ANI_ path']):
    print('Specify the path to the ANI2 potentials in the file, paths.json')
    print('Obtain ANI from its repo')
    print('https://github.com/isayev/ASE_ANI')
    print('Exit Hugin')
    quit()
sys.path.append(paths['ANI_ path'])
from ase_interface import ANIENS
from ase_interface import aniensloader

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

for thefile in filelist:
    outfile = thefile.split('.')[0] + '.xyz'
    if outfile in os.listdir(outpath):
        continue
    ase_atoms = read(inpath + thefile)
    ase_atoms.set_calculator(ANIENS(aniensloader(paths['ANI_ path'] + 'ani_models/ani-2x_8x.info',0,multigpu=False)))
    dyn = LBFGS(ase_atoms)
    dyn.run(fmax=0.001, steps = 2000)
    ase_atoms.write(outpath + outfile, format = 'xyz', plain=True)
