#!/usr/bin/env python

# 2021-04-08

# Kiyoto Aramis Tanemura

# accept the energy file to filter high relative energy structures, and output models enumerated in increasing relative energy.
# Arg1: path to the energy.csv
# Arg2: path to model.smi

import os, sys
import json
import pandas as pd
from ase.io import read 

if len(sys.argv) < 3:
    print('accept the energy file to filter high relative energy structures, and output models enumerated in increasing relative energy.')
    print('Arg1: path to the energy.csv')
    print('Arg2: path to model.smi')
    quit()

with open('arguments.json', 'r') as g:
    args = json.load(g)

adduct = sys.argv[1].split('/')[2]
inpath = '/'.join(sys.argv[1].split('/')[:-1]) + '/'
outpath = '/'.join(sys.argv[1].split('/')[:-2]) + '/charged_model.smi'

# Read the computed single point energy values
df = pd.read_csv(sys.argv[1], index_col = 0)
# retain up to 10 kcal/mol in relative energy
df = df.loc[df['energy'] < 10]
# sort by increasing relative energy
df.sort_values(by='energy', inplace = True)
# retain only up to a given number of models
if df.shape[0] > args['consider_up_to_n_models']:
    df = df.iloc[:args['consider_up_to_n_models']]

with open(sys.argv[2], 'r') as f:
    allsmiles = f.readlines()

for i in [int(x.split('.')[0]) for x in df.index]:
    with open(outpath, 'a') as g:
        g.write(allsmiles[i])
