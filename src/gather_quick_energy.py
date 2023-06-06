#!/usr/bin/env python                                                                                                         

# 2021-04-08

# Kiyoto Aramis Tanemura

# Retrieve single point energy values from QUICK output files
# Arg1: directory to QUICK outputs

import os, sys
import pandas as pd
import numpy as np

if len(sys.argv) < 2:
    print('Retrieve single point energy values from QUICK output files ')
    print('Arg1: directory to QUICK outputs')
    quit()

inpath = sys.argv[1]
if inpath[-1] != '/':
    inpath += '/'

outfiles = [x for x in os.listdir(inpath) if x[-3:] == 'out']

Evals = []

for i in range(len(outfiles)):
    with open(inpath + str(i) + '.out', 'r') as f:
        lines = f.readlines()
    energy = [x for x in lines if 'TOTAL ENERGY' in x]
    energy = energy[0] # get it out of the list
    energy = energy.split()[-1] # it is reported in the last field
    energy = float(energy)
    Evals.append(energy)

min_e = np.min(Evals)
rel_e = [x - min_e for x in Evals]
rel_e = [x * 627.5094740631 for x in rel_e]
df = pd.DataFrame({'energy': rel_e}, index = [str(x) + '.out' for x in range(len(outfiles))])
df.sort_values(by='energy', inplace = True)
df['model'] = range(df.shape[0])
df.to_csv(inpath + 'energy.csv')
print(df)
