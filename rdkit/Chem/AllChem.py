# $Id$
#
#  Copyright (C) 2006-2011  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Import all RDKit chemistry modules

"""
from rdkit import rdBase
from rdkit import RDConfig
from rdkit import DataStructs
from rdkit.Geometry import rdGeometry
from rdkit.Chem import *
from rdkit.Chem.rdPartialCharges import *
from rdkit.Chem.rdDepictor import *
from rdkit.Chem.rdForceFieldHelpers import *
from rdkit.Chem.ChemicalFeatures import *
from rdkit.Chem.rdDistGeom import *
# from rdkit.Chem.rdMolAlign import *
# from rdkit.Chem.MolAlign import *
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdShapeHelpers import *
from rdkit.Chem.rdChemReactions import *
from rdkit.Chem.rdReducedGraphs import *
try:
  from rdkit.Chem.rdSLNParse import *
except ImportError:
  pass
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdqueries import *
from rdkit import ForceField
Mol.Compute2DCoords = Compute2DCoords
Mol.ComputeGasteigerCharges = ComputeGasteigerCharges
import numpy, os
from rdkit.RDLogger import logger
logger = logger()
import warnings


def TransformMol(mol, tform, confId=-1, keepConfs=False):
  """  Applies the transformation (usually a 4x4 double matrix) to a molecule
  if keepConfs is False then all but that conformer are removed

  """
  refConf = mol.GetConformer(confId)
  TransformConformer(refConf, tform)
  if not keepConfs:
    if confId == -1:
      confId = 0
    allConfIds = [c.GetId() for c in mol.GetConformers()]
    for id in allConfIds:
      if not id == confId:
        mol.RemoveConformer(id)
    # reset the conf Id to zero since there is only one conformer left
    mol.GetConformer(confId).SetId(0)


def ComputeMolShape(mol, confId=-1, boxDim=(20, 20, 20), spacing=0.5, **kwargs):
  """ returns a grid representation of the molecule's shape
  """
  res = rdGeometry.UniformGrid3D(boxDim[0], boxDim[1], boxDim[2], spacing=spacing)
  EncodeShape(mol, res, confId, **kwargs)
  return res


def ComputeMolVolume(mol, confId=-1, gridSpacing=0.2, boxMargin=2.0):
  """ Calculates the volume of a particular conformer of a molecule
  based on a grid-encoding of the molecular shape.

  """
  mol = rdchem.Mol(mol)
  conf = mol.GetConformer(confId)
  CanonicalizeConformer(conf)
  box = ComputeConfBox(conf)
  sideLen = (box[1].x - box[0].x + 2 * boxMargin, \
              box[1].y - box[0].y + 2 * boxMargin, \
              box[1].z - box[0].z + 2 * boxMargin)
  shape = rdGeometry.UniformGrid3D(sideLen[0], sideLen[1], sideLen[2], spacing=gridSpacing)
  EncodeShape(mol, shape, confId, ignoreHs=False, vdwScale=1.0)
  voxelVol = gridSpacing ** 3
  occVect = shape.GetOccupancyVect()
  voxels = [1 for x in occVect if x == 3]
  vol = voxelVol * len(voxels)
  return vol


def GenerateDepictionMatching2DStructure(mol, reference, confId=-1, referencePattern=None,
                                         acceptFailure=False, **kwargs):
  """ Generates a depiction for a molecule where a piece of the molecule
     is constrained to have the same coordinates as a reference.

     This is useful for, for example, generating depictions of SAR data
     sets so that the cores of the molecules are all oriented the same
     way.

  Arguments:
    - mol:          the molecule to be aligned, this will come back
                    with a single conformer.
    - reference:    a molecule with the reference atoms to align to;
                    this should have a depiction.
    - confId:       (optional) the id of the reference conformation to use
    - referencePattern:  (optional) an optional molecule to be used to
                         generate the atom mapping between the molecule
                         and the reference.
    - acceptFailure: (optional) if True, standard depictions will be generated
                     for molecules that don't have a substructure match to the
                     reference; if False, a ValueError will be raised

  """
  if reference and referencePattern:
    if reference.GetNumAtoms(onlyExplicit=True) != referencePattern.GetNumAtoms(onlyExplicit=True):
      raise ValueError(
        'When a pattern is provided, it must have the same number of atoms as the reference')
    referenceMatch = reference.GetSubstructMatch(referencePattern)
    if not referenceMatch:
      raise ValueError("Reference does not map to itself")
  else:
    referenceMatch = list(range(reference.GetNumAtoms(onlyExplicit=True)))
  if referencePattern:
    match = mol.GetSubstructMatch(referencePattern)
  else:
    match = mol.GetSubstructMatch(reference)

  if not match:
    if not acceptFailure:
      raise ValueError('Substructure match with reference not found.')
    else:
      coordMap = {}
  else:
    conf = reference.GetConformer()
    coordMap = {}
    for i, idx in enumerate(match):
      pt3 = conf.GetAtomPosition(referenceMatch[i])
      pt2 = rdGeometry.Point2D(pt3.x, pt3.y)
      coordMap[idx] = pt2
  Compute2DCoords(mol, clearConfs=True, coordMap=coordMap, canonOrient=False)


def GenerateDepictionMatching3DStructure(mol, reference, confId=-1, **kwargs):
  """ Generates a depiction for a molecule where a piece of the molecule
     is constrained to have coordinates similar to those of a 3D reference
     structure.

  Arguments:
    - mol:          the molecule to be aligned, this will come back
                    with a single conformer.
    - reference:    a molecule with the reference atoms to align to;
                    this should have a depiction.
    - confId:       (optional) the id of the reference conformation to use

  """
  nAts = mol.GetNumAtoms()
  dm = []
  conf = reference.GetConformer(confId)
  for i in range(nAts):
    pi = conf.GetAtomPosition(i)
    # npi.z=0
    for j in range(i + 1, nAts):
      pj = conf.GetAtomPosition(j)
      # pj.z=0
      dm.append((pi - pj).Length())
  dm = numpy.array(dm)
  Compute2DCoordsMimicDistmat(mol, dm, **kwargs)


def EnumerateLibraryFromReaction(reaction, sidechainSets):
  """ Returns a generator for the virtual library defined by
   a reaction and a sequence of sidechain sets

  >>> from rdkit import Chem
  >>> from rdkit.Chem import AllChem
  >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
  >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
  >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
  >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
  ['CNC=O', 'CCNC=O', 'CNC(C)=O', 'CCNC(C)=O']

  Note that this is all done in a lazy manner, so "infinitely" large libraries can
  be done without worrying about running out of memory. Your patience will run out first:

  Define a set of 10000 amines:
  >>> amines = (Chem.MolFromSmiles('N'+'C'*x) for x in range(10000))

  ... a set of 10000 acids
  >>> acids = (Chem.MolFromSmiles('OC(=O)'+'C'*x) for x in range(10000))

  ... now the virtual library (1e8 compounds in principle):
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[acids,amines])

  ... look at the first 4 compounds:
  >>> [Chem.MolToSmiles(next(r)[0]) for x in range(4)]
  ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']


  """
  if len(sidechainSets) != reaction.GetNumReactantTemplates():
    raise ValueError('%d sidechains provided, %d required' %
                     (len(sidechainSets), reaction.GetNumReactantTemplates()))

  def _combiEnumerator(items, depth=0):
    for item in items[depth]:
      if depth + 1 < len(items):
        v = _combiEnumerator(items, depth + 1)
        for entry in v:
          l = [item]
          l.extend(entry)
          yield l
      else:
        yield [item]

  for chains in _combiEnumerator(sidechainSets):
    prodSets = reaction.RunReactants(chains)
    for prods in prodSets:
      yield prods


def AssignBondOrdersFromTemplate(refmol, mol):
  """ assigns bond orders to a molecule based on the
      bond orders in a template molecule

  Arguments
    - refmol: the template molecule
    - mol: the molecule to assign bond orders to

    An example, start by generating a template from a SMILES
    and read in the PDB structure of the molecule
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("CN1C(=NC(C1=O)(c2ccccc2)c3ccccc3)N")
    >>> mol = AllChem.MolFromPDBFile(os.path.join(RDConfig.RDCodeDir, 'Chem',
    ...          'test_data', '4DJU_lig.pdb'))
    >>> len([1 for b in template.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8
    >>> len([1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    22

    Now assign the bond orders based on the template molecule
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> len([1 for b in newMol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8

    Note that the template molecule should have no explicit hydrogens
    else the algorithm will fail.

    It also works if there are different formal charges (this was github issue 235):
    >>> template=AllChem.MolFromSmiles('CN(C)C(=O)Cc1ccc2c(c1)NC(=O)c3ccc(cc3N2)c4ccc(c(c4)OC)[N+](=O)[O-]')
    >>> mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4FTR_lig.mol'))
    >>> AllChem.MolToSmiles(mol)
    'COC1CC(C2CCC3C(O)NC4CC(CC(O)N(C)C)CCC4NC3C2)CCC1N(O)O'
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> AllChem.MolToSmiles(newMol)
    'COc1cc(-c2ccc3c(c2)Nc2ccc(CC(=O)N(C)C)cc2NC3=O)ccc1[N+](=O)[O-]'

  """
  refmol2 = rdchem.Mol(refmol)
  mol2 = rdchem.Mol(mol)
  # do the molecules match already?
  matching = mol2.GetSubstructMatch(refmol2)
  if not matching:  # no, they don't match
    # check if bonds of mol are SINGLE
    for b in mol2.GetBonds():
      if b.GetBondType() != BondType.SINGLE:
        b.SetBondType(BondType.SINGLE)
        b.SetIsAromatic(False)
    # set the bonds of mol to SINGLE
    for b in refmol2.GetBonds():
      b.SetBondType(BondType.SINGLE)
      b.SetIsAromatic(False)
    # set atom charges to zero;
    for a in refmol2.GetAtoms():
      a.SetFormalCharge(0)
    for a in mol2.GetAtoms():
      a.SetFormalCharge(0)

    matching = mol2.GetSubstructMatches(refmol2, uniquify=False)
    # do the molecules match now?
    if matching:
      if len(matching) > 1:
        logger.warning("More than one matching pattern found - picking one")
      matching = matching[0]
      # apply matching: set bond properties
      for b in refmol.GetBonds():
        atom1 = matching[b.GetBeginAtomIdx()]
        atom2 = matching[b.GetEndAtomIdx()]
        b2 = mol2.GetBondBetweenAtoms(atom1, atom2)
        b2.SetBondType(b.GetBondType())
        b2.SetIsAromatic(b.GetIsAromatic())
      # apply matching: set atom properties
      for a in refmol.GetAtoms():
        a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
        a2.SetHybridization(a.GetHybridization())
        a2.SetIsAromatic(a.GetIsAromatic())
        a2.SetNumExplicitHs(a.GetNumExplicitHs())
        a2.SetFormalCharge(a.GetFormalCharge())
      SanitizeMol(mol2)
      if hasattr(mol2, '__sssAtoms'):
        mol2.__sssAtoms = None  # we don't want all bonds highlighted
    else:
      raise ValueError("No matching found")
  return mol2


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest, sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
