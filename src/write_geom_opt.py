#!/usr/bin/env python

# 2021-04-08

# Kiyoto Aramis Tanemura

# After centroids have been selected by AutoGraph, write out geometry optimization jobs of centroid conformers by QUICK. 
# Arg1: cluster_summary.csv from AutoGraph output
# Arg2: output directory

import os, sys
import pandas as pd
from ase.io import read
from ase.io.xyz import write_xyz

if len(sys.argv) < 3:
    print('After centroids have been selected by AutoGraph, write out geometry optimization jobs of centroid conformers by QUICK.')
    print('Arg1: cluster_summary.csv from AutoGraph output')
    print('Arg2: output directory')
    quit()

inpath = sys.argv[1].split('/')
inpath[3] = 'optimized_conformers'
charge = inpath[2]
inpath = '/'.join(inpath[:-1]) + '/'

outpath = sys.argv[2]
if outpath[-1] != '/':
    outpath += '/'

df = pd.read_csv(sys.argv[1], index_col = 0)
df = df.loc[df['center'] == 1]

for thefile in df.index:
    geom = read(inpath + thefile)
    with open(outpath + thefile.split('.')[0] + '.in', 'w') as f:
        write_xyz(f, [geom])
#    geom.write(outpath + thefile.split('.')[0] + '.in', format='xyz', plain=True)
    with open(outpath + thefile.split('.')[0] + '.in', 'r') as f:
        lines = f.readlines()
    lines[0] = 'DFT B3LYP BASIS=6-31++G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 OPTIMIZE DIPOLE CHARGE=' + charge + '\n'
    lines[1] = '\n'
    with open(outpath + thefile.split('.')[0] + '.in', 'w') as f:
        f.writelines(lines)
