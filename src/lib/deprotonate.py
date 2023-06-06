# 2020-11-15

# Kiyoto Aramis Tanemura

# Deprotonate an ASE Atoms object. Assumes the read molecule is neutral to begin. Must modify if a charged molecules is employed
 
import numpy as np
from ase import Atom
from ase.visualize import view
from ase.geometry.analysis import Analysis
from ase.optimize import LBFGS
from ase.calculators.mopac import MOPAC
from ase.constraints import FixAtoms
from rdkit.Chem.rdmolops import AddHs, GetAdjacencyMatrix

def get_adjacency_matrix(ase_atoms):
    '''Obtain the bond adjacency matrix from ASE Atoms objects. Format to symmetric n x n matrix with no self-loops'''
    ana = Analysis(ase_atoms)
    A = ana.adjacency_matrix[0].A
    A += A.T
    A[range(len(ase_atoms)), range(len(ase_atoms))] = 0
    return A

def deprotonation_sites(ase_atoms, adjacencyMatrix):
    '''Given an ASE structure, return indices of protons which may be abstracted as a list of ints
    Note: we consider only protons bonded to ONFClS, not C regardless of inductive/resonance effects'''
    heteroatoms = [x for x in range(len(ase_atoms)) if ase_atoms.symbols[x] in ['O', 'N', 'F', 'Cl', 'S']]
    deprotonation_sites = [x for x in range(len(ase_atoms)) if ase_atoms.symbols[x] == 'H' and adjacencyMatrix[x,heteroatoms].sum() > 0]
    neighbors = [np.argmax(adjacencyMatrix[x]) for x in deprotonation_sites]
    return deprotonation_sites, neighbors

def deprotonate_structures(ase_atoms, deprotonation_sites, neighbors):
    '''Return list of hypothesized deprotonated structures, as well as the indices where the charge belongs'''
    # Modified 2021-04-07 by KAT. Fix atom coordinates to avoid breaking of compound 
    num_atoms = len(ase_atoms) - 1
    A = get_adjacency_matrix(ase_atoms)
    deprotonation_states = [ase_atoms.copy() for x in range(len(deprotonation_sites))]
    for i in range(len(deprotonation_sites)):
        deprotonation_states[i].pop(deprotonation_sites[i])
        deprotonation_states[i][neighbors[i]].charge = -1
        deprotonation_states[i].set_constraint(FixAtoms(indices=range(num_atoms)))
    # view(deprotonation_states)
    return deprotonation_states

def contains_deprotonation_site(rdkit_mol):
    '''Given a RDKit molecule, return a Boolean whether or not it contains a possible deprotonation site'''
    elements = [x.GetSymbol() for x in AddHs(rdkit_mol).GetAtoms()]
    A = GetAdjacencyMatrix(AddHs(rdkit_mol))
    heteroatoms = [x for x in range(len(elements)) if elements[x] in ['O', 'N', 'F', 'Cl', 'S']]
    H_indices = [x for x in range(len(elements)) if elements[x] == 'H']
    if A[heteroatoms, :][: , H_indices].sum() > 0:
        return True
    else:
        return False
