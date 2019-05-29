import os
import numpy as np
from jpype import java, isJVMStarted, startJVM, getDefaultJVMPath, JPackage
import pycdk

################################## Start JVM #########################################
if not isJVMStarted():
    cdk_path = os.path.join(pycdk.__path__[0], 'cdk-2.2.jar')
    # cdk_path = os.path.join('pycdk/cdk-2.2.jar')
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % cdk_path)
    cdk = JPackage('org').openscience.cdk
else:
    raise OSError ('JVM is already started, please shut down')

############################### Format Conversion #####################################
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

def MolToMOPAC(mol):
    output = java.io.StringWriter()
    writer = cdk.io.program.Mopac7Writer(output)
    writer.write(mol)
    writer.close()
    output = output.toString()
    return output

def MolToSDF(mol):
    output = java.io.StringWriter()
    writer = cdk.io.SDFWriter(output)
    writer.write(mol)
    writer.close()
    output = output.toString()
    return output 

############################# Molecular Properties ##################################
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

############################# Formula and Isotope ##################################
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
    return np.array(output)

def IsotopeSimilarity(isotope_array_1, isotope_array_2, tolerance_ppm=10):
    isotope_1 = IsotopeFromArray(isotope_array_1)
    isotope_2 = IsotopeFromArray(isotope_array_2)
    function = cdk.formula.IsotopePatternSimilarity()
    function.seTolerance(tolerance_ppm)
    output = function.compare(isotope_1, isotope_2)
    return output

def generate_formula(mass, window=0.01, atom_list = {'C': [0, 20], 'H': [0, 20], 'O': [0, 20], 'N': [0, 20], 'P': [0, 20], 'S': [0, 20]}, astring=True):
    ifac = cdk.config.Isotopes.getInstance()
    mfrange = cdk.formula.MolecularFormulaRange()
    builder = cdk.formula.MolecularFormula().getBuilder()
    generator = cdk.formula.MolecularFormulaGenerator
    for atom, (minimum, maximum) in atom_list.items():
        element = ifac.getMajorIsotope(atom)
        mfrange.addIsotope(element, minimum, maximum)
    formula = generator(builder, mass-window, mass+window, mfrange)
    formula = formula.getAllFormulas()
    formula = formula.molecularFormulas()
    if astring:
        formula = [FormulaToString(f) for f in formula]
    return formula

def check_formula(formula, NitrogenRuleCheck=True, RDBERuleCheck=True):
    valid = 1
    if type(formula) == str:
        formula = FormulaFromString(formula)
    if NitrogenRuleCheck:
        checker = cdk.formula.rules.NitrogenRule()
        valid *= checker.validate(formula)
    if RDBERuleCheck:
        checker = cdk.formula.rules.RDBERule()
        valid *= checker.validate(formula)
    if valid > 0:
        return True
    else:
        return False

############################### Fingerprint ########################################
def getFingerprint(mol, fp_type="standard", size=1024, depth=6, transform=True):
    if fp_type == 'maccs':
        nbit = 166
    elif fp_type == 'estate':
        nbit = 79
    elif fp_type == 'pubchem':
        nbit = 881
    elif fp_type == 'klekota-roth':
        nbit = 4860
    else:
        nbit = size
    _fingerprinters = {"standard":cdk.fingerprint.Fingerprinter(size, depth)
                            , "extended":cdk.fingerprint.ExtendedFingerprinter(size, depth)
                            , "substructure": cdk.fingerprint.SubstructureFingerprinter()
                            , "graph":cdk.fingerprint.GraphOnlyFingerprinter(size, depth)
                            , "maccs":cdk.fingerprint.MACCSFingerprinter()
                            , "pubchem":cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
                            , "estate":cdk.fingerprint.EStateFingerprinter()
                            , "hybridization":cdk.fingerprint.HybridizationFingerprinter(size, depth)
                            , "lingo":cdk.fingerprint.LingoFingerprinter(depth)
                            , "klekota-roth":cdk.fingerprint.KlekotaRothFingerprinter()
                            , "shortestpath":cdk.fingerprint.ShortestPathFingerprinter(size)
                            , "signature": cdk.fingerprint.SignatureFingerprinter(depth)
                            , "circular": cdk.fingerprint.CircularFingerprinter()
                            }
    if fp_type in _fingerprinters:
        fingerprinter = _fingerprinters[fp_type]
    else:
        raise IOError('invalid fingerprint type')
    fp = fingerprinter.getBitFingerprint(mol)
    if transform:
        fp = fp.asBitSet()
        bits = []
        idx = fp.nextSetBit(0)
        while idx >= 0:
            bits.append(idx)
            idx = fp.nextSetBit(idx + 1)
        return {'nbit': nbit, 'bits':bits}
    else:
        return fp
    
def TanimotoSimilarity(fingerprint_1, fingerprint_2):
    similarity = cdk.similarity.Tanimoto.calculate(fingerprint_1, fingerprint_2)
    return similarity

################################# Fragmenter #########################################   
def generateFragments(mol, method='MurckoFragmenter', minFragSize=6, singleFrameworkOnly=True, asSmiles=True):
    if method == 'MurckoFragmenter':
        function = cdk.fragment.MurckoFragmenter(singleFrameworkOnly, minFragSize)
    elif method == 'ExhaustiveFragmenter':
        function = cdk.fragment.ExhaustiveFragmenter(minFragSize)
    else:
        raise IOError('Invalid fragmentation method')
    function.generateFragments(mol)
    if asSmiles:
        fragments = function.getFragments()
    else:
        fragments = function.getFragmentsAsContainers()
    return np.array(fragments)

################################# Descriptor #########################################
def getMolecularDescriptorCategories():
    function = cdk.qsar.DescriptorEngine(cdk.qsar.IMolecularDescriptor, cdk.silent.SilentChemObjectBuilder.getInstance())
    return list(function.getAvailableDictionaryClasses())

def getMolecularDescriptor(mol, species='all'):
    function = cdk.qsar.DescriptorEngine(cdk.qsar.IMolecularDescriptor, cdk.silent.SilentChemObjectBuilder.getInstance())
    descriptors = list(function.getDescriptorInstances())
    specifications = list(function.getDescriptorSpecifications())
    categories = []
    for s in specifications:
        try:
            t = list(function.getDictionaryClass(s))
        except:
            t = ['constitutionalDescriptor']
        categories.append(t)
    Descriptors = {}
    if species == 'all':
        species = set(getMolecularDescriptorCategories())
    for i, descriptor in enumerate(descriptors):
        if set(categories[i]).intersection(species) == 0:
            continue
        try:
            name = list(descriptor.getDescriptorNames())[0]
            value = descriptor.calculate(mol).getValue().toString()
            value = value.split(',')
            value = [float(v) for v in value]
        except:
            continue
        Descriptors[name] = value
