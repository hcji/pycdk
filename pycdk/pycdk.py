import os
import numpy as np
from jpype import isJVMStarted, startJVM, getDefaultJVMPath, JPackage
# import pycdk

if not isJVMStarted():
    # cdk_path = os.path.join(pycdk.__path__[0], 'cdk-2.2.jar')
    cdk_path = os.path.join('pycdk/cdk-2.2.jar')
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % cdk_path)
    cdk = JPackage('org').openscience.cdk
else:
    raise OSError ('JVM is already started, please shut down')

def MolFromSmiles(smi):
    function = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
    try:
         mol = function.parseSmiles(smi)
    except:
        raise IOError('invalid smiles input')
    return mol

def MolFromInchi(inchi):
    function = cdk.inchi.InChIGeneratorFactory.getInstance()
    builder = cdk.DefaultChemObjectBuilder.getInstance()
    s = function.getInChIToStructure(inchi, builder)
    mol = s.getAtomContainer()
    return mol

def MolToSmiles(mol):
    function = cdk.smiles.SmilesGenerator(cdk.smiles.SmiFlavor.Isomeric)
    smi = function.create(mol)
    return smi

def MolToInchi(mol):
    function = cdk.inchi.InChIGeneratorFactory.getInstance()
    inchi = function.getInChIGenerator(mol)
    return inchi.getInchi()
    
def MolToInchiKey(mol):
    function = cdk.inchi.InChIGeneratorFactory.getInstance()
    inchi = function.getInChIGenerator(mol)
    return inchi.getInchiKey()

def MolToFormula(mol, string=True):
    function = cdk.tools.manipulator.MolecularFormulaManipulator
    gen = function.getMolecularFormula(mol)
    if string:
        output = function.getString(gen)
    else:
        output = gen
    return output
    
def getMolExactMass(mol):
    function = cdk.tools.manipulator.MolecularFormulaManipulator
    formula = function.getMolecularFormula(mol)
    ExactMass = function.getMajorIsotopeMass(formula)
    return ExactMass

def getMolNaturalMass(mol):
    function = cdk.tools.manipulator.AtomContainerManipulator
    NaturalMass = function.getNaturalExactMass(mol)
    return NaturalMass
    
def getMolTotalFormalCharge(mol):
    function = cdk.tools.manipulator.AtomContainerManipulator
    FormalCharge = function.getTotalFormalCharge(mol)
    return FormalCharge
    
def getMolTotalNegativeFormalCharge(mol):
    function = cdk.tools.manipulator.AtomContainerManipulator
    NegativeFormalCharge = function.getTotalNegativeFormalCharge(mol)
    return NegativeFormalCharge

def getMolTotalPositiveFormalCharge(mol):
    function = cdk.tools.manipulator.AtomContainerManipulator
    PositiveFormalCharge = function.getTotalPositiveFormalCharge(mol)
    return PositiveFormalCharge

def FormulaFromString(string):
    builder = cdk.formula.MolecularFormula().getBuilder()
    formula = cdk.tools.manipulator.MolecularFormulaManipulator.getMolecularFormula(string, builder)
    return formula
    
def FormulaToString(formula):
    string = cdk.tools.manipulator.MolecularFormulaManipulator.getString(formula)
    return string

def getFormulaExactMass(string):
    formula = FormulaFromString(string)
    function = cdk.tools.manipulator.MolecularFormulaManipulator
    ExactMass = function.getMajorIsotopeMass(formula)
    return ExactMass

def getFormulaNaturalMass(string):
    formula = FormulaFromString(string)
    function = cdk.tools.manipulator.MolecularFormulaManipulator
    NaturalMass = function.getNaturalExactMass(formula)
    return NaturalMass

def IsotopeFromString(string, minI=0.01):
    formula = FormulaFromString(string)
    return IsotopeFromFormula(formula, minI)
    
def IsotopeFromFormula(formula, minI=0.01):
    generator = cdk.formula.IsotopePatternGenerator(minI)
    isotopes = generator.getIsotopes(formula)
    isotopes = isotopes.getIsotopes()
    output = [(i.getMass(), i.getIntensity()) for i in isotopes]
    return np.array(output)

def IsotopeFromArray(array):
    isotopes = cdk.formula.IsotopePattern()
    manipulator = cdk.formula.IsotopePatternManipulator
    container = cdk.formula.IsotopeContainer
    for (mass, intensity) in array:
        i = container(mass, intensity)
        isotopes.addIsotope(i)
    output = manipulator.normalize(isotopes)
    output = manipulator.sortByMass(output)
    return output
        
def IsotopeToArray(isotopes):
    isotopes = isotopes.getIsotopes()
    output = [(i.getMass(), i.getIntensity()) for i in isotopes]   
    return output

def IsotopeSimilarity(isotope_array_1, isotope_array_2, tolerance_ppm=10):
    isotope_1 = IsotopeFromArray(isotope_array_1)
    isotope_2 = IsotopeFromArray(isotope_array_2)
    function = cdk.formula.IsotopePatternSimilarity()
    function.seTolerance(tolerance_ppm)
    output = function.compare(isotope_1, isotope_2)
    return output
    
    
