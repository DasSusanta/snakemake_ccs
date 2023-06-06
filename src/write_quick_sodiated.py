#!/usr/bin/env python                                                                                                         

# 2021-12-14

# Kiyoto Aramis Tanemura

# write input to QUICK for monosodiated adducts.
# Arg1: path to output QUICK input script in the sodiated drirectory to be written

import os, sys
import numpy as np
import pandas as pd
from ase.io import read
from ase.io.xyz import write_xyz

if len(sys.argv) < 2:
    print('write input to QUICK for monosodiated adducts.')
    print('Arg1: path to output QUICK input script in the sodiated drirectory to be written')
    quit()

# Get parameters from the path
x, label, x, x, which, model, conf_in = sys.argv[1].split('/')

# get information on the sodiation site
energy = pd.read_csv('results/' + label + '/1/sodiated/site_screen/energy.csv', index_col = 0)
sodiation_index = energy.at[model + '.out', 'sodiation_site'] - 1 # sodiation indices begin at 1, but this is Python

# read input xyz file
geom = read('results/' + label + '/1/ensemble/' + model + '/' + conf_in.replace('in', 'xyz'))
# mutate the proton at the sodiation index to sodium
geom.symbols[sodiation_index] = 'Na'
# get the difference vector of the sodium to center of mass. normalize it to a unit vector
diff_vec = geom[sodiation_index].position - geom.get_center_of_mass()
norm = np.sqrt((diff_vec ** 2).sum())
diff_vec *= 3 / diff_vec
# move the sodium ion by 1 angstrom away from the center of mass
geom[sodiation_index].position += diff_vec

# write out to xyz
with open(sys.argv[1], 'w') as f:
    write_xyz(f, [geom])

# read the xyz and modify the header to QUICK input
with open(sys.argv[1], 'r') as g:
    lines = g.readlines()

if which == 'ensemble':
    lines[0] = 'DFT B3LYP BASIS=6-31++G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 OPTIMIZE DIPOLE CHARGE=1\n'
elif which == 'ensemble_fast':
    lines[0] = 'DFT B3LYP BASIS=6-31++G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 ENERGY DIPOLE CHARGE=1\n'

with open(sys.argv[1], 'w') as h:
    h.writelines(lines)

