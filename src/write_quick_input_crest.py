#!/usr/bin/env python                                                                                                         

# 2021-04-08

# Kiyoto Aramis Tanemura

# write inputs to QUICK for the hypothesized protonation/depotonated molecules by CREST.
# Arg1: directory to input xyz
# Arg2: outdir

import sys
from ase.io import read

if len(sys.argv) < 3:
    print('write inputs to QUICK for the hypothesized protonation/depotonated molecules by CREST.')
    print('Arg1: directory to input xyz')
    print('Arg2: directory to output directory')
    quit()

outpath = sys.argv[2]
if outpath[-1] != '/':
    outpath += '/'

if 'deprotonated' in outpath:
    charge = -1
elif 'protonated' in outpath:
    charge = 1
else:
    print('charge not recognized')
    quit()

geoms = read(sys.argv[1], index=':')
for i in range(len(geoms)):
    geoms[i].write(outpath + str(i) + '.in', format='xyz', plain=True)
    with open(outpath + str(i) + '.in', 'r') as f:
        lines = f.readlines()
    lines[0] = 'DFT B3LYP BASIS=6-311G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 ENERGY CHARGE=' + str(charge) + '\n'
    lines[1] = '\n'
    with open(outpath + str(i) + '.in', 'w') as f:
        f.writelines(lines)
