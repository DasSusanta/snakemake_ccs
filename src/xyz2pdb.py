#!/usr/bin/env python

# 2021-08-09

# Kiyoto Aramis Tanemura

# Take a list of xyz files, convert it to pdb format
# Arg1: Path to clustered_conformers.csv
# Arg2: path to directory of xyz files
# Arg3: path to output directory

import os, sys
import pandas as pd
from ase.io import read

if len(sys.argv) == 1:
    print('After Geometry optimizing structures, retrieve results for subsequent HPCCS calculation. We need energy, xyz coordinates, and Mulliken charges.')
    print('Arg1: Path to clustered_conformers.csv')
    print('Arg2: path to directory of xyz files')
    print('Arg3: path to output directory')
    quit()

inpath = sys.argv[2].split('/')
inpath = '/'.join(inpath[:-1]) + '/'
outpath = sys.argv[3].split('/')
outpath = '/'.join(outpath[:-1]) + '/'

xyz_files = pd.read_csv(sys.argv[1], index_col = 0).index.tolist()
for thefile in xyz_files:
    geom = read(inpath + thefile)
    geom.write(outpath + thefile.replace('xyz', 'pdb'))

