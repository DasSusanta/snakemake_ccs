# 2020-11-15

# Kiyoto Aramis Tanemura

# After obtaining the charged specie in ASE, flattening the molecule back to a mol or SMILES is challenging. Instead, I choose to interface between the same molecule represented in ASE and RDKit to protonate/deprotonate the at the correct site in RDKit. The modified SMILES is then subjected to conformer generation using the ETKDG algorithm implemented in RDKit.
# We featurize the atoms using a simplified implementation of extended connectivity fingerprints (ECFP), used in graph convolutions. Given a bond adjacency matrix and a feature vector (one-hot-encoded elements), we can take the inner product several times to obtain chemical identifiers for each atom. In this manner, we can mapbetween the indices of the molecules in ASE and RDKit.

import numpy as np
from ase.geometry.analysis import Analysis
from rdkit.Chem.rdmolops import GetAdjacencyMatrix, AddHs
from rdkit.Chem import BondType
from .protonate import get_adjacency_matrix

def featurize_ase(ase_atoms):
    '''Return a one-hot-encoder feature matrix to store atom type of each entry'''
    elements = list(ase_atoms.symbols)
    ohe = np.zeros([len(ase_atoms), len(set(ase_atoms.symbols))])
    elements_header = list(set(ase_atoms.symbols))
    for i in range(len(ase_atoms)):
        ohe[i, elements_header.index(elements[i])] = 1
    return ohe, elements_header

def featurize_rdkit(rdkit_mol, elements_header):
    elements = [x.GetSymbol() for x in rdkit_mol.GetAtoms()]
    X = np.zeros([len(elements), len(set(elements))])
    for i in range(len(elements)):
        X[i, elements_header.index(elements[i])] = 1
    return X

def graph_prop(adjacencyMatrix, featureVec, propagation = 5):
    '''Produce ECFP'''
    for i in range(propagation):
        featureVec = adjacencyMatrix @ featureVec
    return featureVec

def find_base_index(protonated_ase_atoms):
    '''Find the index of the base which accepted the proton in protonated ASE Atoms object'''
    A = get_adjacency_matrix(protonated_ase_atoms)
    print(A[-1,:])
    baseindex = np.argmax(A[-1,:])
    return baseindex

    A = GetAdjacencyMatrix(AddHs(mol))
    for i in range(propagation):
        X = A @ X 
    return X

def protonate_X_rdkit(unprotonated_rdkit_mol, rdkit_base_index):
    '''Protonate an heteroatom ONFClS'''
    atom = unprotonated_rdkit_mol.GetAtomWithIdx(rdkit_base_index)
    atom.SetFormalCharge(1)

def protonate_C_rdkit(unprotonated_rdkit_mol, rdkit_base_index):
    '''Protonate an unsaturated C'''
    A = GetAdjacencyMatrix(unprotonated_rdkit_mol)
    neighbors = [x for x in range(A.shape[0]) if A[rdkit_base_index, x] == 1]
    for neighbor in neighbors:
        if unprotonated_rdkit_mol.GetBondBetweenAtoms(rdkit_base_index, neighbor).GetBondTypeAsDouble() > 1:
           unprotonated_rdkit_mol.GetBondBetweenAtoms(rdkit_base_index, neighbor).SetBondType(BondType.SINGLE)
           unprotonated_rdkit_mol.GetAtomWithIdx(neighbor).SetFormalCharge(1)

def protonate_rdkit(unprotonated_ase_atoms, unprotonated_rdkit_mol, baseIndex_ase):
    '''Use a protonated ASE Atoms object as a reference to protonate a molecule in RDKit
    RDKit molecule should have explicit H specified'''
    A_ase = get_adjacency_matrix(unprotonated_ase_atoms)
    X_ase, elements_header = featurize_ase(unprotonated_ase_atoms)

    A_rdkit = GetAdjacencyMatrix(unprotonated_rdkit_mol)
    X_rdkit = featurize_rdkit(unprotonated_rdkit_mol, elements_header)

    X_ase = graph_prop(A_ase, X_ase)
    X_rdkit = graph_prop(A_rdkit, X_rdkit)

    reference = X_ase[baseIndex_ase, :]
    comparison = X_rdkit - reference
    rdkit_base_index = np.sum(comparison, axis = 1).tolist().index(0.0)

    # Modified 2020-11-17 by KAT. Treat protonation of heteroatom or carbon differently
    elements = [x.GetSymbol() for x in unprotonated_rdkit_mol.GetAtoms()]
    if elements[rdkit_base_index] in ['O', 'N', 'F', 'Cl', 'S']:
        protonate_X_rdkit(unprotonated_rdkit_mol, rdkit_base_index)
    elif elements[rdkit_base_index] == 'C':
        protonate_C_rdkit(unprotonated_rdkit_mol, rdkit_base_index)
#    atom = unprotonated_rdkit_mol.GetAtomWithIdx(rdkit_base_index)
#    atom.SetFormalCharge(1)

    return unprotonated_rdkit_mol

def deprotonate_rdkit(protonated_ase_atoms, protonated_rdkit_mol, acidIndex_ase):
    '''Use information from the deprotonated ASE Atoms object to place appropriate negative
    charge on the RDKit molecule'''
    A_ase = get_adjacency_matrix(protonated_ase_atoms)
    X_ase, elements_header = featurize_ase(protonated_ase_atoms)

    A_rdkit = GetAdjacencyMatrix(protonated_rdkit_mol)
    X_rdkit = featurize_rdkit(protonated_rdkit_mol, elements_header)
    
    X_ase = graph_prop(A_ase, X_ase)
    X_rdkit = graph_prop(A_rdkit, X_rdkit)
 
    reference = X_ase[acidIndex_ase]
    comparison = X_rdkit - reference
    rdkit_acid_index = np.sum(comparison, axis = 1).tolist().index(0.0)
    neighbor_indices = [x for x in range(A_rdkit.shape[0]) if A_rdkit[rdkit_acid_index, x] == 1]
    elements = [x.GetSymbol() for x in protonated_rdkit_mol.GetAtoms()]
    proton_index = [x for x in neighbor_indices if elements[x] == 'H']

    atom = protonated_rdkit_mol.GetAtomWithIdx(rdkit_acid_index)
    atom.SetFormalCharge(-1)

    return protonated_rdkit_mol

















