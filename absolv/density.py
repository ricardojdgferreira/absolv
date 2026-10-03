#!/usr/bin/python

# every argument of this script should describe a molecule
# either as a MDL .mol file or a SMILES string.
# Please note that the RDKit library is required.

import sys, os
from collections import Counter
import rdkit.Chem.rdMolDescriptors
from rdkit import Chem

# volume increments for the model
Vinc = {
    'H1-0': 0,      # H volume dumped into heavy atoms
    'Al6-3': 7.27,  # This is how RDkit sees Al3-0
    'B3-0': 7.61,
    'B3-1': 22.09,
    'Br1-0': 27.52,
    'C1-0': 14.68,
    'C2-0': 9.89,
    'C2-1': 22.77,
    'C3-0': -0.00,
    'C3-1': 13.23,
    'C3-2': 27.17,
    'C4-0': -8.40,
    'C4-1': 4.46,
    'C4-2': 16.57,
    'C4-3': 29.58,
    'Cl1-0': 24.74,
    'F1-0': 17.63,
    'Ge4-0': 9.36,
    'Ge4-1': 25.92,
    'I1-0': 35.64,
    'N1-0': 14.09,
    'N2-0': 7.42,
    'N2-1': 18.14,
    'N3-0': -3.08,
    'N3-1': 7.74,
    'N3-2': 17.81,
    'O1-0': 14.89,
    'O2-0': 6.25,
    'O2-1': 11.78,
    'P3-0': 10.42,
    'P4-0': -1.94,
    'P4-1': 10.06,
    'Pb4-0': 14.29,
    'R5': 9.41,
    'R6': 6.89,
    'R<5': 10.89,
    'R>6': 3.75,
    'S1-0': 25.92,
    'S2-0': 14.90,
    'S2-1': 26.14,
    'S3-0': 5.58,
    'S4-0': -3.74,
    'Se2-0': 19.00,
    'Si4-0': 9.28,
    'Si4-1': 23.35,
    'Sn4-0': 14.47,
    'Ti4-0': 6.09,
    'aromat': 1.82
}

# Set of parameters fitted using <10 entries 
rares = set(('Al6-3', 'B3-1', 'C1-0', 'Ge4-1', 'N2-1', 'P4-1', 'Pb4-0', 'S3-0', 'Se2-0', 'Ti4-0'))

def GetARH(atom):
    """ return an atom code consistent with the keys of the 'Vinc' dictionary """ 
    # be careful to take into account nD = number of deuterium atoms !
    nD = len([x for x in atom.GetNeighbors() if x.GetMass()>2 and x.GetSymbol()=='H'])
    if nD > 0:
        print('Warning: nD %u %s' %(nD, arg))
    return atom.GetSymbol() + str(atom.GetTotalDegree()) + '-' + str(atom.GetTotalNumHs()+nD)

def areAromatic(atoms):
    """ test whether all input atoms (typically from a ring) are aromatic """
    for at in atoms:
        if not at.GetIsAromatic():
            return False
    return True

def setRingData(dic, mol, maxi=6, mini=5):
    """ add keys relative to rings in dictionary 'dic' for molecule 'mol' """
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        size = len(ring)
        if size > maxi:
            label = 'R>'+str(maxi)
        elif size < mini:
            label = 'R<'+str(mini)
        else:
            label = 'R%u' %len(ring)
        dic[label] += 1
        # contribute also +1 aromatic ring ?
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if areAromatic(atoms):
            dic['aromat'] += 1
    return dic

def GetMol(arg):
    """ return RDkit molecule from input argument: either .mol file or SMILES """
    if not os.path.isfile(arg):
        try:
            return Chem.MolFromSmiles(arg)
        except ValueError:
            return None
    if arg.endswith('.sdf'):
        return Chem.MolFromMolFile(arg)
    return None

def Density(arg):            # iterate over all input arguments
    mol = GetMol(arg)               # each one is for a molecule
    if mol is None:
        print('Error: compound %s cannot be read' %arg)
    mf = Counter([a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms()])  # empirical formula
    if mf['C'] == 0:
        print('Warning: %s not in applicability domain (no C atom)' %arg)
    if mf['H'] < 2:
        print('Warning: %s exhibits %u H atom (low accuracy expected)' %(arg, mf['H']))
    dic = Counter()                 # counts number of fragments/rings
    for at in mol.GetAtoms():
        k = GetARH(at)              # 'k' defines the kind of fragment associated with each non-H atom
        dic[k] +=1                  # the number of such fragments is incremented
    setRingData(dic, mol)           # add the number of aromatic rings to 'dic'
    dangerous = set(dic.keys()) & rares  # do we use parameters adjusted using <10 entries ?
    if dangerous:
        print('Warning: possibly ill-defined parameters for %s: ' %arg + ' '.join(sorted(list(dangerous))))
    try:
        Vm = sum([dic[k]*Vinc[k] for k in dic])  # adds increments involved
    except KeyError:                # 'k' not found in Vinc for unknown fragments
        print(dic)
        print('Error: missing parameter(s) for compound: ' +arg)        # nevertheless print a line in output file
    if Vm == 0:                     # compound made only of H/D atoms
        print('Error: cannot calculate molar volume for compound %s' %arg)
    Mw = rdkit.Chem.rdMolDescriptors.CalcExactMolWt(mol)  # molecular weight needed...
    calcDensity = Mw/Vm                                   # to convert molar volume to density
    # return only density
    return float(f"{calcDensity:.4f}")
