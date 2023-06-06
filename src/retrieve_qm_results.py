#!/usr/bin/env python

# 2021-04-10

# Kiyoto Aramis Tanemura

# After Geometry optimizing structures, retrieve results for subsequent HPCCS calculation. We need energy, xyz coordinates, and Mulliken charges.
# Arg1: Path to the QUICK output file

import os, sys
import pandas as pd

if len(sys.argv) == 1:
    print('After Geometry optimizing structures, retrieve results for subsequent HPCCS calculation. We need energy, xyz coordinates, and Mulliken charges.')
    print('Arg1: Path to the QUICK output file')
    quit()

inpath = sys.argv[1].split('/')
inpath = '/'.join(inpath[:-1]) + '/'

energy_keyword = ' TOTAL ENERGY '
geometry_keyword = 'GEOMETRY INFORMATION'
charge_keyword = 'MULLIKEN'
atomic_radii = {'C': '1.70', 'H': '1.20', 'O': '1.52', 'N': '1.55', 'F': '1.47', 'S': '1.80', 'Cl': '1.75', 'NA': '2.27', 'P': '1.80'}

# read QUICK output file
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
energy_index = [x for x in range(len(lines)) if energy_keyword in lines[x]][-1] # get last total energy from the cylce of geom opt.
e_val = float(lines[energy_index].split()[-1]) * 627.5094740631

# select index of first occurance of geometry keyword.
geom_index = [x for x in range(len(lines)) if geometry_keyword in lines[x]][0]
# Drop all above the first xyz entry, which occurs four lines below keyword.
lines = lines[geom_index + 4:]
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

# geometric information is recorded. Next, record Mulliken charges
charge_index = [x for x in range(len(lines)) if charge_keyword in lines[x]][0]
charge_list = []
for i in range(len(elems)):
    charge_list.append(lines[charge_index + 1 + i].split()[1])

xs = [x.rjust(7) for x in xs]
ys = [x.rjust(7) for x in ys]
zs = [x.rjust(7) for x in zs]
charge_list = [x.rjust(7) for x in charge_list]
# We have retrived the relevant information from file: energy, coordinates, and Mulliken charge
# Produce pqr file from coordinates and charge. Write out energy in its own csv file.
atom = ['ATOM' for x in range(len(elems))]
serial = [str(x) for x in range(1, len(elems) + 1)]
residue = ['MOL' for x in range(len(elems))]
chain = ['1' for x in range(len(elems))]
radii = [atomic_radii[x] for x in elems]

pqr = pd.DataFrame({'atom': atom, 'serial': serial, 'element': elems, 'residue': residue, 'chain': chain, 'X': xs, 'Y': ys, 'Z': zs, 'charge': charge_list, 'rradius': radii})
pqr.to_csv(sys.argv[1][:-3] + 'pqr', header = False, index = False, sep = '\t')

xyz = [str(elems[i]) + '\t' + str(xs[i]) + '\t' + str(ys[i]) + '\t' + str(zs[i]) + '\n' for i in range(len(xs))]
with open(sys.argv[1][:-3] + 'xyz', 'w') as g:
    g.writelines([str(len(xs)) + '\n', '\n'] + xyz)

with open(inpath + 'energy.csv', 'a') as h:
    h.write(sys.argv[1].split('/')[-1] + '\t' + str(e_val) + '\n')
