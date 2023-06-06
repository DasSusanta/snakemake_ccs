#!/usr/bin/env python
# 2021-12-14

# Kiyoto Aramis Tanemura

# write QUICK input file for single point calculation for monosodiated adduct. 
# Arg1: path to QUICK input file to be written in the sodiated/site_screen directory

import sys
import pandas as pd

if len(sys.argv) == 1:
    print('write QUICK input file for single point calculation for monosodiated adduct.')
    print('Arg1: path to QUICK input file to be written in the sodiated/site_screen directory')
    quit()

# read label and model from outpath
outpath = sys.argv[1]
x, label, x, x, x, model_in = outpath.split('/')
# obtain the sodiation site as the proton with the greatest charge
path_to_directory_of_protonation_adduct_ensemble = '/'.join(outpath.split('/')[:3]) + '/ensemble/' + model_in.split('.')[0] + '/'
considered_confoermers_df = pd.read_csv(path_to_directory_of_protonation_adduct_ensemble + 'energy_considered_conformers.csv', index_col = 0)
min_e_conformer_name = considered_confoermers_df.at[0, 'file'].replace('out', 'pqr')
min_e_conformer = pd.read_csv(path_to_directory_of_protonation_adduct_ensemble + min_e_conformer_name, delimiter='\t', names = ['atom', 'num', 'elem', 'mol', 'chain', 'x', 'y', 'z', 'charge', 'radius'])
sodiation_site = min_e_conformer.loc[min_e_conformer['elem'] == 'H']['charge'].idxmax()
sodiation_site = min_e_conformer.at[sodiation_site, 'num'] - 1 #subtract 1 since index begins at 1, and we want to work with Python index system

# open the xyz of the min_e_conformer
with open(path_to_directory_of_protonation_adduct_ensemble + min_e_conformer_name.replace('pqr', 'xyz'), 'r') as f:
    lines = f.readlines()

# sodiate by replacing the proton at the sodiation_index to a sodium atom. skip first two lines of the header
lines[sodiation_site + 2] = lines[sodiation_site + 2].replace('H', 'Na')

# enter arguments to QUICK. We only consider monoprotonated species for now. 
lines[0] = 'DFT B3LYP BASIS=6-31++G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 ENERGY CHARGE=1\n'

with open(sys.argv[1], 'w') as g:
    g.writelines(lines)
