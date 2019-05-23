# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:03:22 2019

@author: hcji
"""

# from pycdk import *

# read mol from inchi/smiles
smi = 'c1ccc(cc1)CN(c2cc(ccc2[N+](=O)[O-])c3c(nc(nc3CC)N)N)C'
inchi = 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
mol = MolFromSmiles(smi)
mol = MolFromInchi(inchi)

# get inchi/smiles
smi = MolToSmiles(mol)
inchi = MolToInchi(mol)
inchikey = MolToInchiKey(mol)

# generate sdf/mopac
mopac = MolToMOPAC(mol)
sdf = MolToSDF(mol)

# get informations from mol
MolToFormula(mol)
getMolExactMass(mol)
getMolNaturalMass(mol)
getMolTotalFormalCharge(mol)
getMolTotalNegativeFormalCharge(mol)
getMolTotalPositiveFormalCharge(mol)

# get fingerprint and compare similarity
mol1 = MolFromSmiles('CCCO')
mol2 = MolFromSmiles('COCC')
fingerprint = getFingerprint(mol1, fp_type="standard", size=1024, depth=6, transform=True)
fingerprint1 = getFingerprint(mol1, fp_type="standard", size=1024, depth=6, transform=False)
fingerprint2 = getFingerprint(mol2, fp_type="standard", size=1024, depth=6, transform=False)
TanimotoSimilarity(fingerprint1, fingerprint2)

# formula from/to string
string = 'C2H5OH'
formula = FormulaFromString(string)
string1 = FormulaToString(formula)

# get mass from formula string
getFormulaExactMass(string)
getFormulaNaturalMass(string)

# isotope from string/formula
isotope_array = IsotopeFromString(string, minI=0.01)
isotope_array = IsotopeFromFormula(formula, minI=0.01)

# cdk.isotope from array
isotope = IsotopeFromArray(isotope_array)

# generate formula from mass
mass = 18.03383
window = 0.01
atom_list = {'C': [0, 20], 'H': [0, 20], 'O': [0, 20], 'P': [0, 20], 'S': [0, 20]}
formulas = generate_formula(mass, window, atom_list)

# check if formula is valid
valid = check_formula('C5H16O')

# compare isotope pattern
string1 = 'C5H16O'
string2 = 'C4H12O2'
isotope_array_1 = IsotopeFromString(string1, minI=0.001)
isotope_array_2 = IsotopeFromString(string2, minI=0.001)
similarity1 = IsotopeSimilarity(isotope_array_1, isotope_array_2, 10)
similarity2 = IsotopeSimilarity(isotope_array_1, isotope_array_2, 3.4)