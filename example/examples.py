# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:03:22 2019

@author: hcji
"""

# from pycdk import *

# read mol from inchi/smiles
smi = 'CCCO'
inchi = 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
mol = MolFromSmiles(smi)
mol = MolFromInchi(inchi)

# get inchi/smiles
smi = MolToSmiles(mol)
inchi = MolToInchi(mol)
inchikey = MolToInchiKey(mol)

# get informations from mol
MolToFormula(mol)
getMolExactMass(mol)
getMolNaturalMass(mol)
getMolTotalFormalCharge(mol)
getMolTotalNegativeFormalCharge(mol)
getMolTotalPositiveFormalCharge(mol)

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

# compare isotope pattern
string1 = 'C5H16O'
string2 = 'C4H12O2'
isotope_array_1 = IsotopeFromString(string1, minI=0.001)
isotope_array_2 = IsotopeFromString(string2, minI=0.001)
similarity1 = IsotopeSimilarity(isotope_array_1, isotope_array_2, 10)
similarity2 = IsotopeSimilarity(isotope_array_1, isotope_array_2, 3.4)