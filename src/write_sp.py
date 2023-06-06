#!/usr/bin/env python

# 2021-08-11

# Kiyoto Aramis Tanemura

# Write single point calculation inputs for conformers selected as centroids
# Arg1: path to directory containing centroid conformer files
# Arg2: path to directory to save input files

import os, sys
from ase.io import read

if len(sys.argv) < 3:
    print('Write single point calculation inputs for conformers selected as centroids')
    print('Arg1: path to directory containing centroid conformer files')
    print('Arg2: path to directory to save input files')
    quit()

inpath = '/'.join(sys.argv[1].split('/')) + '/'
outpath = '/'.join(sys.argv[2].split('/')) + '/'

adduct = inpath.split('/')[2]

filelist = [x for x in os.listdir(inpath) if x[-3:] == 'xyz']

for thefile in filelist:
    geom = read(inpath + thefile)
    geom.write(outpath + thefile.replace('xyz', 'in'), format = 'xyz')
    with open(outpath + thefile.replace('xyz', 'in'), 'r') as f:
        lines = f.readlines()
    lines[0] = 'DFT B3LYP BASIS=6-311G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 ENERGY DIPOLE CHARGE=' + adduct + '\n'
    lines[1] = '\n'
    with open(outpath + thefile.replace('xyz', 'in'), 'w') as g:
        g.writelines(lines)


