import os
import jpype
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

def MolToSmiles(mol):
    function = cdk.smiles.SmilesGenerator(cdk.smiles.SmiFlavor.Isomeric)
    smi = function.create(mol)
    return smi
    
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
    

