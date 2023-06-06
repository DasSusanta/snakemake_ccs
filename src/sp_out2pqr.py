#!/usr/bin/env python                                                                                                         

# 2021-08-26

# Kiyoto Aramis Tanemura

# Take single point output file. Convert it to pqr file.
# Arg1: single point output file

import os, sys
import pandas as pd
import numpy as np

if len(sys.argv) < 2:
    print('Take single point output file. Convert it to pqr file.')
    print('Arg1: single point output file')
    quit()

thefile = sys.argv[1]

energy_keyword = 'TOTAL ENERGY'

geometry_keyword = '-- INPUT GEOMETRY --'
charge_keyword = 'MULLIKEN'
atomic_radii = {'C': '1.70', 'H': '1.20', 'O': '1.52', 'N': '1.55', 'F': '1.47', 'S': '1.80', 'Cl': '1.75', 'NA': '2.27'}

with open(thefile, 'r') as f:
    lines = f.readlines()

# obtain single point energy value
energy_line = [x for x in lines if energy_keyword in x]
energy_line = energy_line[0]
energy_val = float(energy_line.split('=')[1]) * 627.5094740631 # kcal/mol

# select index of first occurance of geometry keyword.  
geom_index = [x for x in range(len(lines)) if geometry_keyword in lines[x]][0]
# Drop all above the first xyz entry
lines = lines[geom_index+1:]
i = 0
elems = []
xs = []
ys = []
zs = []
while len(lines[i].split()) == 4:
    element, X, Y, Z = lines[i].split()
    elems.append(element)
    xs.append(X)
    ys.append(Y)
    zs.append(Z)
    i += 1

with open(thefile, 'r') as g:
    lines = g.readlines()
charge_index = [x for x in range(len(lines)) if charge_keyword in lines[x]][0]
charge_list = []
for i in range(len(elems)):
    charge_list.append(lines[charge_index + 1 + i].split()[1])

xs = [x.rjust(7) for x in xs]
ys = [x.rjust(7) for x in ys]
zs = [x.rjust(7) for x in zs]
charge_list = [x.rjust(7) for x in charge_list]

# Produce pqr file from coordinates and charge.
atom = ['ATOM' for x in range(len(elems))]
serial = [str(x) for x in range(1, len(elems) + 1)]
residue = ['MOL' for x in range(len(elems))]
chain = ['1' for x in range(len(elems))]
radii = [atomic_radii[x] for x in elems]

# write out PQR file
pqr = pd.DataFrame({'atom': atom, 'serial': serial, 'element': elems, 'residue': residue, 'chain': chain, 'X': xs, 'Y': ys, 'Z': zs, 'charge': charge_list, 'radius': radii})
pqr.to_csv(thefile.replace('out', 'pqr'), header = False, index = False, sep = '\t')

# Lettus also write out the XYZ files
xyz_lines = [str(len(elems)) + '\n', '\n'] + ['\t'.join([elems[i], xs[i], ys[i], zs[i]]) + '\n' for i in range(len(elems))]
with open(thefile.replace('out', 'xyz'), 'w') as h:
    h.writelines(xyz_lines)

with open('/'.join(sys.argv[1].split('/')[:-1]) + '/energy.csv', 'a') as i:
    i.write(sys.argv[1].split('/')[-1] + '\t' + str(energy_val) + '\n')
