#!/usr/bin/env python
# 2021-12-14

# Kiyoto Aramis Tanemura

# Single point energy values are collected for sodiation isomers. Gather QUICK output files for energy and summarize in one file named energy.csv
# Arg1: path to output files in sodiated directory

import os, sys
import pandas as pd
import numpy as np

inpath = sys.argv[1]
if inpath[-1] != '/':
    inpath += '/'

outfiles = [x for x in os.listdir(inpath) if x[-3:] == 'out'] # model0.out, model1.out, etc

Evals = []
sodiation_indices = []

for thefile in outfiles:
    with open(inpath + thefile, 'r') as f:
        lines = f.readlines()
    energy = [x for x in lines if 'TOTAL ENERGY' in x]
    energy = energy[0] # get it out of the list
    if '=' in energy: # I think I included this to remove the "=" if the digits are large and touch the equals sign
        energy = energy.replace('=', '')
    energy = energy.split()[-1] # it is reported in the last field  
    energy = float(energy)
    Evals.append(energy)

    # let us also record the sodiation site.
    path_to_directory_of_protonation_adduct_ensemble = '/'.join(inpath.split('/')[:3]) + '/ensemble/' + thefile.split('.')[0] + '/'
    considered_confoermers_df = pd.read_csv(path_to_directory_of_protonation_adduct_ensemble + 'energy_considered_conformers.csv', index_col = 0)
    min_e_conformer_name = considered_confoermers_df.at[0, 'file'].replace('out', 'pqr')
    min_e_conformer = pd.read_csv(path_to_directory_of_protonation_adduct_ensemble + min_e_conformer_name, delimiter='\t', names = ['atom', 'num', 'elem', 'mol', 'chain', 'x', 'y', 'z', 'charge', 'radius'])
    sodiation_site = min_e_conformer.loc[min_e_conformer['elem'] == 'H']['charge'].idxmax()
    sodiation_site = min_e_conformer.at[sodiation_site, 'num']
    sodiation_indices.append(sodiation_site)

rel_e = [x - np.min(Evals) for x in Evals]
rel_e = [x * 627.5094740631 for x in rel_e]
df = pd.DataFrame({'energy': [x * 627.5094740631 for x in Evals], 'rel_energy': rel_e, 'sodiation_site': sodiation_indices}, index = [x.split()[0] for x in outfiles])
df.sort_values(by='energy', inplace = True)
df.to_csv(inpath + 'energy.csv')
print(df)
