#!/usr/bin/env python                                                                                                                                                                                

# 2021-04-06

# Kiyoto Aramis Tanemura

# Once conformers are generated, refine them using ANI-2x potential.
# Arg1: input mol file
# Arg2: output xyz file

import os, sys
from ase.io import read
from lib.ANI2_refinement import ani_minimize

if len(sys.argv) < 3:
    print('Refine conformer file by ANI-2x potential')
    print('Arg1: input mol file')
    print('Arg2: output xyz file')
    quit()
conf = read(sys.argv[1])
ani_minimize(conf)
conf.write(sys.argv[2], format = 'xyz', plain=True)
