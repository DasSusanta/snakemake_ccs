# 2020-11-19

# Kiyoto Aramis Tanemura

# Automate conformational ensemble generation for predicted protonated/deprotonated states of small orgnaic molecules

import os, sys
import json
import pandas as pd
import numpy as np
# Atomic Simulation Environment
from ase.io import read
# RDKit
from rdkit.Chem import MolToMolBlock, SmilesMolSupplier
from rdkit.Chem.rdmolops import GetFormalCharge, AddHs, RemoveHs
from rdkit.Chem.rdmolfiles import MolFromMolFile, MolFromSmiles, MolToSmiles

with open('paths.json', 'r') as f:
    paths = json.load(f)

if not os.path.isfile(paths['xyz2mol_path']):
    print('Specify path to xyz2mol in paths.json')
    print('Follow download/installation from its repository:')
    print('https://github.com/jensengroup/xyz2mol')
    print('Exit Hugin program')
    quit()

if not os.path.isdir(paths['AutoGraph_path']):
    print('Specify the path to AutoGraph in paths.json')
    print('Contact Kiyoto to obtain the package')
    print('tanemur1@msu.edu')
    print('Exit Hugin program')
    quit()

sys.path.append(paths['AutoGraph_path'])
from AutoGraph import AutoGraph

# Custom

from lib.cleanSMILES import neutralize_sanitize
from lib.confGen import generateconformations, getMMFFenergiesConf, saveConfMol
from lib.protonate import protonate_structures, optimize_structures, retrieve_min_E_state, contains_protonation_site
from lib.deprotonate import get_adjacency_matrix, deprotonation_sites, deprotonate_structures, contains_deprotonation_site
from lib.ase_rdkit_interface import find_base_index, protonate_rdkit, deprotonate_rdkit
from lib.ANI2_refinement import ani_minimize

class Hugin:
    '''Workflow for ensemble generation of molecular ions to be used in CCS prediction, 
    eventually to be completed with the ab initio CCS prediction step'''
    def __init__(self, outpath = 'Hugin', label = 'Hugin'):
        '''outpath: (string) path to output directory. Can be made already or added
        label: (string) identifier for the particular job
        If not supplied, find based on structural features.'''
        self.outpath = outpath
        if self.outpath[-1] !='/':
            self.outpath += '/'
        self.label = label
        if not os.path.isdir(self.outpath):
            os.mkdir(self.outpath)

    def check_adducts(self, rdk_mol):
        '''Decide applicable adducts. This should return candidate adducts to query'''
        charge_initial = GetFormalCharge(rdk_mol)
        if charge_initial == 1:
            return [0]
        elif charge_initial == 0:
            bronsted_base = contains_protonation_site(rdk_mol)
            bronsted_acid = contains_deprotonation_site(rdk_mol)
            if bronsted_base:
                if bronsted_acid:
                    return [1, -1]
                return [1]
            elif bronsted_acid:
                return[-1]
        else:
            print('The supplied molecule has a charge that cannot yet be handled by this class. \
                  Please check the molecule if it was interpreted correctly:', rdk_mol)

    def protonation_protocol(self, ase_atoms, max_models):
        '''Protonate an ASE Atoms object and find likely protonation state. Protonate RDKit mol at corresponding position'''
        prot_states = protonate_structures(ase_atoms)
        list_ase_atoms, energies = optimize_structures(prot_states, 1)
        os.system('touch ' + self.outpath + 'prot_model.smi')
        for i in range(np.min([max_models, len(energies)])):
            list_ase_atoms[energies.index(np.sort(energies)[i])].write(self.outpath + 'prot_model' + str(i) + '.xyz', plain = True)
            os.system('python ' + paths['xyz2mol_path'] + ' ' + self.outpath + 'prot_model' + str(i) + '.xyz --charge 1 >> ' + self.outpath + 'prot_model.smi')
        suppl = SmilesMolSupplier(self.outpath + 'prot_model.smi',delimiter='\n',titleLine=False)
        return suppl

    def deprotonation_protocol(self, ase_atoms, max_models):
        '''Deprotonate an ASE Atoms object and find likely protonation state. 
        Derotonate RDKit mol at corresponding position'''
        A = get_adjacency_matrix(ase_atoms)
        sites, neighbor = deprotonation_sites(ase_atoms, A)
        deprot_states = deprotonate_structures(ase_atoms, sites, neighbor)
        list_ase_atoms, energies = optimize_structures(deprot_states, -1)
        os.system('touch ' + self.outpath + 'dep_model.smi')
        for i in range(np.min([max_models, len(energies)])):
            list_ase_atoms[energies.index(np.sort(energies)[i])].write(self.outpath + 'dep_model' + str(i) + '.xyz', plain= True)
            os.system('python ' + paths['xyz2mol_path'] + ' ' + self.outpath + 'dep_model' + str(i) + '.xyz --charge -1 >> ' + self.outpath + 'dep_model.smi')
        suppl = SmilesMolSupplier(self.outpath + 'dep_model.smi',delimiter='\n',titleLine=False)
        return suppl

    def run(self, smiles, adducts = [], conformers = 1000, charge_models=1):
        '''Perform full workflow, from smiles to ensemble structure
        smiles: (string) SMILES of the molecule
        adducts: (int or list of ints) The particular adducts to query from [1, -1, 0] for ["[M+H]+", "[M-H]-", "[M]+"] respectively
        conformers: (int) maximum number of conformers to generate
        prot_models: (int) max number of predicted protonation model to consider'''
        print('---Begin Hugin Program---')
        rdk_mol = MolFromSmiles(smiles)
        rdk_mol = neutralize_sanitize(rdk_mol)
        if type(adducts) == type(0):
            adducts = [adducts] # if an integer is passed, place in list
        elif len(adducts) == 0: # If adducts are not supplied, use function to determine the adducts to query
            adducts = self.check_adducts(rdk_mol)
        # Generate ten conformers, and choose minimum energy conformer. Save as mol
        print('Read the SMILES,', smiles)
        print('Generating one representative conformer')
        confs = generateconformations(rdk_mol, addH = True)
        energies = getMMFFenergiesConf(confs)
        confids = [x.GetId() for x in confs.GetConformers()]
        minEconf = confids[np.argmin(energies)]
        print(MolToMolBlock(confs, confId=minEconf),file=open(self.outpath + self.label + '_starting_conformer.mol','w+'))
        ase_atoms = read(self.outpath + self.label + '_starting_conformer.mol')
        print('Conformer was generated and read. Now modeling possible charge states')
        ionized_species = [[],[],[]] # record RDKit molecules of the respective ionized species.
        if 1 in adducts:
            ionized_species[1] = self.protonation_protocol(ase_atoms, charge_models)
        if -1 in adducts:
            ionized_species[-1] = self.deprotonation_protocol(ase_atoms, charge_models)
        if 0 in adducts:
            ionized_species[0].append(ase_atoms)
        dirnames = ['molecular_ion', 'protonated', 'deprotonated']
        print('Charge states have been determined. Now generate conformers for each.')
        for i in adducts:
            if not os.path.isdir(self.outpath + dirnames[i]):
                os.mkdir(self.outpath + dirnames[i])
            print('---', dirnames[i], '---')
            # Iterate through retained predicted ionized species
            for j in range(len(ionized_species[i])):
                print('---model', str(j), '---')
                print('Generating ', conformers, ' conformers')
                if not os.path.isdir(self.outpath + dirnames[i] + '/model' + str(j)):
                    os.mkdir(self.outpath + dirnames[i] + '/model' + str(j))
                confs = generateconformations(ionized_species[i][j], numConf = conformers, addH = True)
                saveConfMol(confs, self.outpath + dirnames[i] + '/model' + str(j) + '/generated_conformers/')

                if not os.path.isdir(self.outpath + dirnames[i] + '/model' + str(j) + '/optimized_conformers/'):
                    os.mkdir(self.outpath + dirnames[i] + '/model' + str(j) + '/optimized_conformers/')
                conf_files = [x for x in os.listdir(self.outpath + dirnames[i] + '/model' + str(j) + '/generated_conformers/') if x[-3:] == 'mol']

                print('Refining ', len(conf_files), ' conformers using ANI2-x potentials')
                confs = [read(self.outpath + dirnames[i] + '/model' + str(j) + '/generated_conformers/' + x) for x in conf_files]
                conf_files = [x.replace('mol', 'xyz') for x in conf_files]
                for k in range(len(conf_files)):
                    ani_minimize(confs[k])
                    confs[k].write(self.outpath + dirnames[i] + '/model' + str(j) + '/optimized_conformers/' + conf_files[k], plain=True)
                # The corresponding ANI energy values may be important. Record them as energy.csv
                energies = [x.get_potential_energy() for x in confs]
                pd.DataFrame({'energy': energies}, index = conf_files).to_csv(self.outpath + dirnames[i] + '/model' + str(j) + '/optimized_conformers/' + 'energy.csv')
        
                # Cluster conformers using AutoGraph, selecting centroids by minimum energy structures
                if not os.path.isdir(self.outpath + dirnames[i] + '/model' + str(j) + '/clustered_conformers/'):
                    os.mkdir(self.outpath + dirnames[i] + '/model' + str(j) + '/clustered_conformers/')
                ag = AutoGraph()
                ag.run(self.outpath + dirnames[i] + '/model' + str(j) + '/optimized_conformers/', self.outpath + dirnames[i] + '/model' + str(j) + '/clustered_conformers/', self.outpath + dirnames[i] + '/model' + str(j) + '/optimized_conformers/energy.csv')
