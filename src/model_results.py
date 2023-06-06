#!/usr/bin/env python

# 2021-04-10

# Kiyoto Aramis Tanemura

# After HPCCS results have been obtained, gather them for the model.
# Arg1: path to the HPCCS outputs.

import os, sys
import pandas as pd
import numpy as np
from ase.io import read
from ase.build import minimize_rotation_and_translation

if len(sys.argv) == 1:
    print('After HPCCS results have been obtained, gather them for the model.')
    print('Arg1: path to the HPCCS outputs.')
    quit()

kB = 0.001985875 # Boltzmann constant in kcal/(mol K)
T = 298 # room temperature in Kelvin                
beta = 1 / (kB * T)
def boltzmann_contribution(rel_E_values):
    '''Provided a list of relative energy values (kcal/mol), return the Boltzmann weighted average energy'''
    Bfactors = [np.exp(- beta * x) for x in rel_E_values]
    Z = np.sum(Bfactors)
    return [x / Z for x in Bfactors]

wdpath = sys.argv[1]
if wdpath[-1] == '/':
    wdpath += '/'

energy = pd.read_csv(wdpath + 'energy.csv', delimiter = '\t', names = ['file', 'energy'])
#energy['conf_id'] = [x.split('.')[0] for x in energy['file']]
#energy.set_index('conf_id', drop = True, inplace = True)
energy.dropna(inplace = True)
energy['energy'] = pd.to_numeric(energy['energy'])
energy['rel_energy'] = energy['energy'] - energy['energy'].min()
energy.sort_values(by = 'rel_energy', inplace = True)
energy.reset_index(inplace = True, drop = True)

# Modified 2021-05-08 by KAT. I noticed some sampled conformers converged to nearly identical conformers upon geometry optimization
# The redundancy in sampling can introduce erroneous result when calculating the ensemble average. 
# Let us remove any higher energy structure redundant within 0.1 Å to any lower or equal energy structure.
for i in range(energy.shape[0] - 1, -1, -1): # start from highest energy structure
    sys0 = read(wdpath + energy.at[i, 'file'].replace('out', 'xyz'))
    for j in range(i - 1, -1, -1):
        sys1 = read(wdpath + energy.at[j, 'file'].replace('out', 'xyz'))
        minimize_rotation_and_translation(sys0, sys1)
        diff = sys0.positions - sys1.positions
        rmsd = np.sqrt(np.mean((diff ** 2).sum(axis = 1)))
        if rmsd <= 0.1:
            energy.drop(i, inplace = True)
            print('dropped')
            break

energy.to_csv(wdpath + 'energy_considered_conformers.csv')

ccs_list = []
std_list = []

for thefile in [x.replace('out', 'hpccs') for x in energy['file']]:
    with open(wdpath + thefile, 'r') as f:
        lines = f.readlines()
    ccsline = [x for x in lines if 'Average TM cross section' in x][0]
    stdline = [x for x in lines if 'Standard deviation (percent)' in x][0]
    ccs_list.append(float(ccsline.split()[-1]))
    std_list.append(float(stdline.split()[-1]))

energy['ccs'] = ccs_list
energy['rel_std'] = [x / 100 for x in std_list] # get relative standard deviation as fraction [0,1]
energy['fraction'] = boltzmann_contribution(energy['rel_energy'])

final_ccs = (energy['ccs'] * energy['fraction']).sum() / energy['fraction'].sum()
final_abs_std = np.sqrt((energy['rel_std'] * energy['ccs'] * energy['fraction']).pow(2).sum())
final_rel_std = final_abs_std/final_ccs

with open(wdpath + 'ccs.txt', 'w') as g:
    g.write('ccs:\t' + str(final_ccs) + '\nabsolute error:\t' + str(final_abs_std)+ '\nrelative error:\t' + str(final_rel_std))

energy.to_csv(wdpath + 'results.csv')
print(energy)
print('ccs:\t' + str(final_ccs) + '\nabsolute error:\t' + str(final_abs_std)+ '\nrelative error:\t' + str(final_rel_std))
