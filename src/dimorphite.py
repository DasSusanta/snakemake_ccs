#!/usr/bin/env python

# 2021-04-20

# Kiyoto Aramis Tanemura

# It will be advantageous to use an empirical pKa prediction algorithm to produce charged models in CCS prediction. dimorphite-dl returns the expected model for a given pH range, but does not exactly perform the protonation site screening. This code attempts to use dimorphite-dl, but to output all likely protonated/deprotonated states.

import os, sys
import json
from rdkit import Chem
from rdkit.Chem.rdmolops import GetFormalCharge
import pandas as pd 

with open('paths.json', 'r') as f:
    paths = json.load(f)

with open('arguments.json', 'r') as g:
    args = json.load(g)

sys.path.append(paths['dimorphite'])
import dimorphite_dl

# read input SMILES
df = pd.read_csv('data/input.smi', delimiter = '\t', names = ['hmdb_id', 'smiles'])
#df['molecule'] = [Chem.MolFromSmiles(s) for s in df['smiles']]

# for each of the molecules,
for i in range(df.shape[0]):
    if not os.path.isdir('results/' + df.at[i, 'hmdb_id']):
        os.mkdir('results/' + df.at[i, 'hmdb_id'])
    mol = Chem.MolFromSmiles(df.at[i, 'smiles'])
    site_screen = dimorphite_dl.run_with_mol_list([mol], min_ph=args['pH_range'][0], max_ph = args['pH_range'][1], max_variants = args['max_variants'])
    charges = [GetFormalCharge(m) for m in site_screen]
    charge_set = set(charges)
    for theCharge in charge_set:
        if theCharge not in args['consider_adducts_of_charges']:
            continue
        if not os.path.isdir('results/' + df.at[i, 'hmdb_id'] + '/' + str(theCharge)):
            os.mkdir('results/' + df.at[i, 'hmdb_id'] + '/' + str(theCharge))
        if not os.path.isfile('results/' + df.at[i, 'hmdb_id'] + '/' + str(theCharge) + '/model.smi'):
            charge_index = [x for x in range(len(charges)) if charges[x] == theCharge]
            with open('results/' + df.at[i, 'hmdb_id'] + '/' + str(theCharge) + '/model.smi', 'w') as f:
                f.writelines([Chem.MolToSmiles(m) + '\n' for m in [site_screen[x] for x in charge_index]])
