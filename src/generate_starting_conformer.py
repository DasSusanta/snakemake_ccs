#!/usr/bin/env python3

# 2021-04-03

# Kiyoto Aramis Tanemura

# Generate starting conformer for protonation state prediction. Run as executable.
# Arg1: input smiles 
# Arg2: output path

import os, sys
import numpy as np
from rdkit.Chem import MolToMolBlock#, SmilesMolSupplier
from rdkit.Chem.rdmolfiles import MolFromSmiles
from lib.cleanSMILES import neutralize_sanitize
from lib.confGen import generateconformations, getMMFFenergiesConf, saveConfMol

if len(sys.argv) == 1:
    print('Prepare starting geometry to input to Hugin')
    print('arguments (in order): SMILES output_path')
    quit()

smiles = sys.argv[1]
outpath = '' # if output path isn't provided, output in current directory.

if len(sys.argv) > 2:
    outpath = sys.argv[2]
    if outpath[-1] != '/':
        outpath += '/'

if not os.path.isdir(outpath):
    os.mkdir(outpath)

rdk_mol = MolFromSmiles(smiles)
rdk_mol = neutralize_sanitize(rdk_mol)
confs = generateconformations(rdk_mol, addH = True)
energies = getMMFFenergiesConf(confs)
confids = [x.GetId() for x in confs.GetConformers()]
minEconf = confids[np.argmin(energies)]
print(MolToMolBlock(confs, confId=minEconf),file=open(outpath + 'starting_conformer.mol','w+'))
with open(outpath + 'smiles.smi', 'w') as f:
    f.write(smiles)
