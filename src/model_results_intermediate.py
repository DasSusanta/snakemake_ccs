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
if wdpath[-1] != '/':
    wdpath += '/'

conf_path = wdpath.replace('ensemble_fast', 'clustered_conformers') + 'centroids/'
filelist = [x for x in os.listdir(wdpath) if 'hpccs'in x]

energy = pd.read_csv(wdpath + 'energy.csv', delimiter = '\t', names = ['file', 'energy'])
energy.dropna(inplace = True)
energy['energy'] = pd.to_numeric(energy['energy'])
energy['relE'] = energy['energy'] - energy['energy'].min()
energy = energy.loc[energy['relE'] < 5].copy()
energy.sort_values(by = 'relE', inplace = True)
energy.reset_index(inplace = True, drop = True)

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
energy['fraction'] = boltzmann_contribution(energy['relE'])

final_ccs = (energy['ccs'] * energy['fraction']).sum() / energy['fraction'].sum()
final_abs_std = np.sqrt((energy['rel_std'] * energy['ccs'] * energy['fraction']).pow(2).sum())
final_rel_std = final_abs_std/final_ccs

with open(wdpath + 'ccs.txt', 'w') as g:
    g.write('ccs:\t' + str(final_ccs) + '\nabsolute error:\t' + str(final_abs_std)+ '\nrelative error:\t' + str(final_rel_std))

energy.to_csv(wdpath + 'results.csv')
print(energy)
print('ccs:\t' + str(final_ccs) + '\nabsolute error:\t' + str(final_abs_std)+ '\nrelative error:\t' + str(final_rel_std))
