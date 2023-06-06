#!/usr/bin/env python

# 2021-04-05

# Kiyoto Aramis Tanemura

# Given a starting mol file, find protonation/deprotonation sites and produce models for the ion. The script will look for smiles.smi and starting_conformer.mol in the working directory path
# Arg1: the working directory path to obtain smiles.smi and starting_conformer.mol

import os, sys
import json
import numpy as np
from ase.io import read
from rdkit.Chem.rdmolfiles import MolFromMolFile, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import GetFormalCharge, AddHs, RemoveHs

from lib.protonate import protonate_structures, optimize_structures, retrieve_min_E_state, contains_protonation_site
from lib.deprotonate import get_adjacency_matrix, deprotonation_sites, deprotonate_structures, contains_deprotonation_site

with open('paths.json', 'r') as f:
    paths = json.load(f)

if len(sys.argv) == 1:
    print('This executable script produces ionized (protonated/deprotonated) model of the starting geometry')
    print('Arg1: the working directory path to obtain smiles.smi and starting_conformer.mol')
    quit()

wdpath = sys.argv[1]
if wdpath[-1] != '/':
    wdpath += '/'

# read the starting geometry and SMILES
ase_atoms = read(wdpath + 'starting_conformer.mol')
with open(wdpath + 'smiles.smi', 'r') as f:
    smiles = f.read()
smiles = smiles.strip()
rdk_mol = MolFromSmiles(smiles)

# Obtain possible models (molecular ion (m), protonated(p), or deprotonated(d))
adducts = []
charge_initial = GetFormalCharge(rdk_mol)
if charge_initial == 1:
    adducts.append('m')
elif charge_initial == 0:
    if contains_protonation_site(rdk_mol):
        adducts.append('p')
    if contains_deprotonation_site(rdk_mol):
        adducts.append('d')

# Perform the protonation protocol
if 'p' in adducts:
    protpath = wdpath + 'protonated/'
    if not os.path.isdir(protpath):
        os.mkdir(protpath)
    prot_states = protonate_structures(ase_atoms) # list of protonated structures
    list_ase_atoms, energies = optimize_structures(prot_states, 1, label = protpath) # optimize structures with charge of 1
    os.system('touch ' + protpath + 'model.smi')
    rank = np.argsort(energies).tolist() # get the ranking of energies
    if not os.path.isdir(protpath + 'model_xyz/'):
        os.mkdir(protpath + 'model_xyz/')
    for i in range(len(energies)):
        # write out protonated models as xyz from low to high energy
        list_ase_atoms[rank.index(i)].write(protpath + 'model_xyz/model' + str(i) + '.xyz', plain = True)
        # append the smiles for corresponding structure to prot_model.smi
        os.system(paths['xyz2mol_path'] + ' ' + protpath + 'model_xyz/model' + str(i) + '.xyz --charge 1 >> ' + protpath + 'model.smi')    

if 'd' in adducts:
    deprpath = wdpath + 'deprotonated/'
    if not os.path.isdir(deprpath):
        os.mkdir(deprpath)
    A = get_adjacency_matrix(ase_atoms)
    sites, neighbor = deprotonation_sites(ase_atoms, A)
    deprot_states = deprotonate_structures(ase_atoms, sites, neighbor)
    list_ase_atoms, energies = optimize_structures(deprot_states, -1, label = deprpath)
    os.system('touch ' + deprpath + 'model.smi')
    rank = np.argsort(energies).tolist() # get the ranking of energies
    if not os.path.isdir(deprpath + 'model_xyz/'):
        os.mkdir(deprpath + 'model_xyz/')
    for i in range(len(energies)):
        # write out deprotonated models as xyz from low to high energy
        list_ase_atoms[rank.index(i)].write(deprpath + 'model_xyz/model' + str(i) + '.xyz', plain = True)
        # append the smiles for corresponding structure to model.smi                                                        
        os.system(paths['xyz2mol_path'] + ' ' + deprpath + 'model_xyz/model' + str(i) + '.xyz --charge -1 >> ' + deprpath + 'model.smi')

if 'm' in adducts:
    if not os.path.isdir(wdpath + 'molecular_ion/'):
        os.mkdir(wdpath + 'molecular_ion/')
    if not os.path.isdir(wdpath + 'molecular_ion/model_xyz/'):
        os.mkdir(wdpath + 'molecular_ion/model_xyz/')
    with open(wdpath + 'molecular_ion/model.smi', 'w') as f:
        f.write(smiles)
    ase_atoms.write(wdpath + 'molecular_ion/model_xyz/model0.xyz', plain = True)
