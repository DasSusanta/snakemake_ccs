#!/usr/bin/env python

# 2021-04-06

# Kiyoto Aramis Tanemura

# After refining geometries with ANI-2x, calculate the single point energy for each.
# Arg1: path to directory of ANI optimized structures.

import os, sys
import json
import pandas as pd
from ase.io import read
import torchani

with open('paths.json', 'r') as f:
    paths = json.load(f)

if len(sys.argv) < 2:
    print('After refining geometries with ANI-2x, calculate the single point energy for each.')
    print('Arg1: path to directory of ANI optimized structures.')
    quit()

path = sys.argv[1]
if path[-1] != '/':
    path += '/'

#calculator = torchani.models.ANI2x().ase()
##SDAS
calculator = torchani.models.ANI2x(periodic_table_index=True).double()


E_list = []
conf_files = [x for x in os.listdir(path) if x[-3:] == 'xyz'] # read only xyz files
for thefile in conf_files:
    geom = read(path + thefile)
    #geom.set_calculator(calculator)
    geom.set_calculator(calculator.ase())
    E_list.append(geom.get_potential_energy())

pd.DataFrame({'energy': E_list}, index = conf_files).to_csv(path + 'energy.csv')
