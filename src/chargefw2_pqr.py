#!/usr/bin/env python

# 2021-08-09

# Kiyoto Aramis Tanemura

# Take the empirical charge calculation output by ChargeFW2, and write out PQR file.
# Arg1: path to directory with relevant files

import os, sys
import pandas as pd

if len(sys.argv) == 1:
    print('Take the empirical charge calculation output by ChargeFW2, and write out PQR file.')
    print('Arg1: path to directory with relevant files')
    quit()

inpath = sys.argv[1].split('/')
inpath = '/'.join(inpath[:-1]) + '/'

pdbfiles = [x for x in os.listdir(inpath) if x[-3:] == 'pdb']
charges = [x for x in os.listdir(inpath) if x[-3:] == 'txt']

atomic_radii = {'C': '1.70', 'H': '1.20', 'O': '1.52', 'N': '1.55', 'F': '1.47', 'S': '1.80', 'Cl': '1.75'}

for thefile in pdbfiles:
    with open(inpath + thefile, 'r') as f:
        lines = f.readlines()
    lines = [x.split() for x in lines if 'ATOM' in x]
    elems = [x[2] for x in lines]
    xs = [x[5] for x in lines]
    ys = [x[6] for x in lines]
    zs = [x[7] for x in lines]

    with open(inpath + thefile + '.txt', 'r') as f:
        cont = f.readlines()
    charge_list = cont[-1].split()

    xs = [x.rjust(7) for x in xs]
    ys = [x.rjust(7) for x in ys]
    zs = [x.rjust(7) for x in zs]
    charge_list = [x.rjust(7) for x in charge_list]
    # Produce pqr file from coordinates and charge. Write out energy in its own csv file.
    atom = ['ATOM' for x in range(len(elems))]
    serial = [str(x) for x in range(1, len(elems) + 1)]
    residue = ['MOL' for x in range(len(elems))]
    chain = ['1' for x in range(len(elems))]
    radii = [atomic_radii[x] for x in elems]

    pqr = pd.DataFrame({'atom': atom, 'serial': serial, 'element': elems, 'residue': residue, 'chain': chain, 'X': xs, 'Y': ys, 'Z': zs, 'charge': charge_list, 'radius': radii})
    pqr.to_csv(inpath + thefile.replace('pdb', 'pqr'), header = False, index = False, sep = '\t')
