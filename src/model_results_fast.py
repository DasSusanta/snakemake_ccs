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

kB = 8.617333262145 * 10 ** (-5) # Boltzmann constant in eV / K
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
energy = pd.read_csv(wdpath + 'energy.csv', index_col = 0)
energy['rel_energy'] = energy['energy'] - energy['energy'].min()
energy.sort_values(by = 'rel_energy', inplace = True)

i = 0
while i < energy.shape[0]:
    sys0 = read(wdpath + energy.iloc[i].name.replace('xyz', 'pdb'))
    drop_indices = []
    for j in range(i + 1, energy.shape[0]):
        sys1 = read(wdpath + energy.iloc[j].name.replace('xyz', 'pdb'))
        minimize_rotation_and_translation(sys0, sys1)
        diff = sys0.positions - sys1.positions
        rmsd = np.sqrt(np.mean((diff ** 2).sum(axis = 1)))
        if rmsd <= 0.1:
            drop_indices.append(j)
    energy.drop(energy.iloc[drop_indices].index, inplace = True)
    i += 1

ccs_list = []
std_list = []

for thefile in [x.replace('xyz', 'hpccs') for x in energy.index]:
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