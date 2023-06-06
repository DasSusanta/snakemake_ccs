#!/usr/bin/env python

# 2021-04-06

# Kiyoto Aramis Tanemura

# Cluster ANI-2x optimized conformers using AutoGraph.
# Arg1: path to directory containing ANI optimized structures
# Arg2: path of the output directory
# Arg3: path to the energy.csv

import os, sys
import json
with open('paths.json', 'r') as f:
    paths = json.load(f)
sys.path.append(paths['AutoGraph_path'])
from AutoGraph import AutoGraph

if len(sys.argv) < 4:
    print('Cluster ANI-2x optimized conformers using AutoGraph.')
    print('Arg1: path to directory containing ANI optimized structures')
    print('Arg2: path of the output directory')
    print('Arg3: path to the energy.csv')
    quit()

ag = AutoGraph(copy_conformers = False, randomize = True)
ag.run(sys.argv[1], sys.argv[2], sys.argv[3])
