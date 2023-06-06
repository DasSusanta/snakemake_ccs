# 2020-11-15

# Kiyoto Aramis Tanemura

# Refine generated files using ANI2 potentials

import os, sys
import json
from ase.io import read
from ase.optimize import LBFGS
from ase.md.langevin import Langevin
from ase import units
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

def ani_minimize(ase_atoms, fmax=0.001, steps = 2000):
    '''Use ANI2-x potentials to geometry optimize an ASE Atoms object'''
    ase_atoms.set_calculator(ANIENS(aniensloader(paths['ANI_ path'] + 'ani_models/ani-2x_8x.info',0,multigpu=False)))
    dyn = LBFGS(ase_atoms)
    dyn.run(fmax=fmax, steps = steps)

