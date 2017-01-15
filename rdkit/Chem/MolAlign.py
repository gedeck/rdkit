'''
Created on 12 Jan 2017

@author: peter
'''
import warnings

import numpy

from rdkit.Chem.ChemicalForceFields import UFFGetMoleculeForceField
from rdkit.Chem.rdDistGeom import EmbedMolecule

# atom based alignments
from rdkit.Chem.rdMolAlign import AlignmentParameters, AlignMol
from rdkit.Chem.rdMolAlign import AlignMolConformers, GetAlignmentTransform
from rdkit.Chem.rdMolAlign import GetBestRMS, GetConformerRMS
# O3A alignment
from rdkit.Chem.rdMolAlign import GetCrippenO3A, GetCrippenO3AForProbeConfs
from rdkit.Chem.rdMolAlign import GetO3A, GetO3AForProbeConfs
from rdkit.Chem.rdMolAlign import O3A, RandomTransform


def pyGetBestRMS(ref, probe, refConfId=-1, probeConfId=-1, maps=None):
  """ Returns the optimal RMS for aligning two molecules, taking
  symmetry into account. As a side-effect, the probe molecule is
  left in the aligned state.

  Arguments:
    - ref: the reference molecule
    - probe: the molecule to be aligned to the reference
    - refConfId: (optional) reference conformation to use
    - probeConfId: (optional) probe conformation to use
    - maps: (optional) a list of lists of (probeAtomId,refAtomId)
      tuples with the atom-atom mappings of the two molecules.
      If not provided, these will be generated using a substructure
      search.

  Note:
  This function will attempt to align all permutations of matching atom
  orders in both molecules, for some molecules it will lead to 'combinatorial
  explosion' especially if hydrogens are present.
  Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing the
  atom order.

  """
  if not maps:
    matches = ref.GetSubstructMatches(probe, uniquify=False)
    if not matches:
      raise ValueError('mol %s does not match mol %s' % (ref.GetProp('_Name'),
                                                         probe.GetProp('_Name')))
    if len(matches) > 1e6:
      warnings.warn("{} matches detected for molecule {}, this may lead to a performance slowdown.".
                    format(len(matches), probe.GetProp('_Name')))
    maps = [list(enumerate(match)) for match in matches]
  bestRMS = 1000.
  for amap in maps:
    rms = AlignMol(probe, ref, prbCid=probeConfId, refCid=refConfId, atomMap=amap)
    if rms < bestRMS:
      bestRMS = rms
      bestMap = amap

  # finally repeat the best alignment :
  if bestMap != amap:
    AlignMol(probe, ref, prbCid=probeConfId, refCid=refConfId, atomMap=bestMap)
  return bestRMS


def pyGetConformerRMS(mol, confId1, confId2, atomIds=None, prealigned=False):
  """ Returns the RMS between two conformations.
  By default, the conformers will be aligned to the first conformer
  of the molecule (i.e. the reference) before RMS calculation and,
  as a side-effect, will be left in the aligned state.

  Arguments:
    - mol:        the molecule
    - confId1:    the id of the first conformer
    - confId2:    the id of the second conformer
    - atomIds:    (optional) list of atom ids to use a points for
                  alingment - defaults to all atoms
    - prealigned: (optional) by default the conformers are assumed
                  be unaligned and will therefore be aligned to the
                  first conformer

  """
  # align the conformers if necessary
  # Note: the reference conformer is always the first one
  if not prealigned:
    if atomIds:
      AlignMolConformers(mol, confIds=[confId1, confId2], atomIds=atomIds)
    else:
      AlignMolConformers(mol, confIds=[confId1, confId2])

  # calculate the RMS between the two conformations
  conf1 = mol.GetConformer(id=confId1)
  conf2 = mol.GetConformer(id=confId2)
  ssr = 0
  for i in range(mol.GetNumAtoms()):
    d = conf1.GetAtomPosition(i).Distance(conf2.GetAtomPosition(i))
    ssr += d * d
  ssr /= mol.GetNumAtoms()
  return numpy.sqrt(ssr)


def GetConformerRMSMatrix(mol, atomIds=None, prealigned=False):
  """ Returns the RMS matrix of the conformers of a molecule.
  As a side-effect, the conformers will be aligned to the first
  conformer (i.e. the reference) and will left in the aligned state.

  Arguments:
    - mol:     the molecule
    - atomIds: (optional) list of atom ids to use a points for
               alignment - defaults to all atoms
    - prealigned: (optional) by default the conformers are assumed
                  be unaligned and will therefore be aligned to the
                  first conformer

  Note that the returned RMS matrix is symmetrically, i.e. it is the
  lower half of the matrix, e.g. for 5 conformers:
  rmsmatrix = [ a,
                b, c,
                d, e, f,
                g, h, i, j]
  This way it can be directly used as distance matrix in e.g. Butina
  clustering.

  """
  # if necessary, align the conformers
  # Note: the reference conformer is always the first one
  rmsvals = []
  if not prealigned:
    if atomIds:
      AlignMolConformers(mol, atomIds=atomIds, RMSlist=rmsvals)
    else:
      AlignMolConformers(mol, RMSlist=rmsvals)
  else:  # already prealigned
    for i in range(1, mol.GetNumConformers()):
      rmsvals.append(GetConformerRMS(mol, 0, i, atomIds=atomIds, prealigned=prealigned))
  # loop over the conformations (except the reference one)
  cmat = []
  for i in range(1, mol.GetNumConformers()):
    cmat.append(rmsvals[i - 1])
    for j in range(1, i):
      cmat.append(GetConformerRMS(mol, i, j, atomIds=atomIds, prealigned=True))
  return cmat


def ConstrainedEmbed(mol, core, useTethers=True, coreConfId=-1, randomseed=2342,
                     getForceField=UFFGetMoleculeForceField, **kwargs):
  """ generates an embedding of a molecule where part of the molecule
  is constrained to have particular coordinates

  Arguments
    - mol: the molecule to embed
    - core: the molecule to use as a source of constraints
    - useTethers: (optional) if True, the final conformation will be
        optimized subject to a series of extra forces that pull the
        matching atoms to the positions of the core atoms. Otherwise
        simple distance constraints based on the core atoms will be
        used in the optimization.
    - coreConfId: (optional) id of the core conformation to use
    - randomSeed: (optional) seed for the random number generator


    An example, start by generating a template with a 3D structure:
    >>> from rdkit.Chem import AllChem
    >>> from rdkit.Chem import MolAlign
    >>> template = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1")
    >>> AllChem.EmbedMolecule(template)
    0
    >>> AllChem.UFFOptimizeMolecule(template)
    0

    Here's a molecule:
    >>> mol = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")

    Now do the constrained embedding
    >>> newmol = MolAlign.ConstrainedEmbed(mol, template)

    Demonstrate that the positions are the same:
    >>> newp=newmol.GetConformer().GetAtomPosition(0)
    >>> molp=mol.GetConformer().GetAtomPosition(0)
    >>> list(newp-molp)==[0.0,0.0,0.0]
    True
    >>> newp=newmol.GetConformer().GetAtomPosition(1)
    >>> molp=mol.GetConformer().GetAtomPosition(1)
    >>> list(newp-molp)==[0.0,0.0,0.0]
    True

  """
  match = mol.GetSubstructMatch(core)
  if not match:
    raise ValueError("molecule doesn't match the core")
  coordMap = {}
  coreConf = core.GetConformer(coreConfId)
  for i, idxI in enumerate(match):
    corePtI = coreConf.GetAtomPosition(i)
    coordMap[idxI] = corePtI

  ci = EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, **kwargs)
  if ci < 0:
    raise ValueError('Could not embed molecule.')

  algMap = [(j, i) for i, j in enumerate(match)]

  if not useTethers:
    # clean up the conformation
    ff = getForceField(mol, confId=0)
    for i, idxI in enumerate(match):
      for j in range(i + 1, len(match)):
        idxJ = match[j]
        d = coordMap[idxI].Distance(coordMap[idxJ])
        ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.)
    ff.Initialize()
    n = 4
    more = ff.Minimize()
    while more and n:
      more = ff.Minimize()
      n -= 1
    # rotate the embedded conformation onto the core:
    rms = AlignMol(mol, core, atomMap=algMap)
  else:
    # rotate the embedded conformation onto the core:
    rms = AlignMol(mol, core, atomMap=algMap)
    ff = getForceField(mol, confId=0)
    conf = core.GetConformer()
    for i in range(core.GetNumAtoms()):
      p = conf.GetAtomPosition(i)
      pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
      ff.AddDistanceConstraint(pIdx, match[i], 0, 0, 100.)
    ff.Initialize()
    n = 4
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
      more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
      n -= 1
    # realign
    rms = AlignMol(mol, core, atomMap=algMap)
  mol.SetProp('EmbedRMS', str(rms))
  return mol


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import sys
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
