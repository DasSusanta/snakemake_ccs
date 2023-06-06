# 2020-11-15

# Kiyoto Aramis Tanemura

# Take input molecule as SMILES. Use RDKit functions to neturalize and standardize molecules.

from rdkit import Chem

def neutralize_atoms(mol): # obtained from RDKit Cookbook
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def neutralize_sanitize(mol):
    mol_out = neutralize_atoms(mol)
    Chem.SanitizeMol(mol_out)
    return mol_out

