#!/usr/bin/env python                                                                                                         

# 2021-04-08

# Kiyoto Aramis Tanemura

# write inputs to QUICK for the hypothesized protonation/depotonated molecules by CREST.
# Arg1: directory to input xyz
# Arg2: outdir

import os, sys
from ase.io import read

if len(sys.argv) < 3:
    print('write inputs to QUICK for the hypothesized protonation/depotonated molecules by PM3/COSMO.')
    print('Arg1: input directory')
    print('Arg2: output directory')
    quit()

inpath = sys.argv[1]
if inpath[-1] != '/':
    inpath += '/'

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

infiles = [x for x in os.listdir(inpath) if x[-3:] == 'xyz']
for thefile in infiles:
    geom = read(inpath + thefile)
    geom.write(outpath + thefile.split('.')[0] + '.in', format='xyz', plain=True)
    with open(outpath + thefile.split('.')[0] + '.in', 'r') as f:
        lines = f.readlines()
    lines[0] = 'DFT B3LYP BASIS=6-311G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 ENERGY CHARGE=' + str(charge) + '\n'
    lines[1] = '\n'
    with open(outpath + thefile.split('.')[0] + '.in', 'w') as f:
        f.writelines(lines)
