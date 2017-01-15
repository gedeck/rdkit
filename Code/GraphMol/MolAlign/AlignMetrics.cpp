//
//  Copyright (C) 2001-2016 Peter Gedeck
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AlignMolecules.h"
#include <math.h>
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

/*
 * since there’s no “ignore hs while matching” argument, I would copy both
 * molecules, remove the Hs, find the best RMS, and then use the transform from
 * that alignment to transform the original molecule
 */

// template <class T, class U>
// std::vector<std::vector<U> > ConvertVectVect(T vectvect) {
//  std::vector<std::vector<U> > vect;
//  size_t size = python::extract<unsigned int>(vectvect.attr("__len__")());
//  vect.resize(size);
//  for (size_t i = 0; i < size; ++i) {
//    unsigned int len1 =
//        python::extract<unsigned int>(vectvect[i].attr("__len__")());
//    std::vector<U> &v = vect[i];
//    v.reserve(len1);
//    for (unsigned int j = 0; j < len1; ++j) {
//      U u = python::extract<U>(vectvect[i][j]);
//      v.push_back(u);
//    }
//  }
//  return vect;
//}
// std::vector<MatchVectType> atomMaps = ConvertVectVec<python::list,
// MatchType>(maps);

namespace RDKit {
namespace MolAlign {

double getBestRMS2(const ROMol &refMol, ROMol &prbMol, int refCid, int prbCid,
                   bool includeHydrogens,
                   std::vector<MatchVectType> *atomMaps) {
  RDGeom::Point3DConstPtrVect refPoints, prbPoints;

  if (!atomMaps || atomMaps->size() == 0) {
    // we have to figure out all mappings between the two molecule
    const bool uniquify = true;
    const bool recursionPossible = true;
    const bool useChirality = false;
    const bool useQueryQueryMatches = true;
    unsigned int nmatches =
        SubstructMatch(refMol, prbMol, *atomMaps, uniquify, recursionPossible,
                       useChirality, useQueryQueryMatches);
    if (nmatches == 0) {
      throw MolAlignException(
          "No sub-structure match found between the probe and query mol");
    }
  }

  double bestRMS = 10000.0;
  MatchVectType bestMap;
  bool bestAtLast = false;
  BOOST_FOREACH (const MatchVectType &atomMap, *atomMaps) {
    double rms = alignMol(prbMol, refMol, prbCid, refCid, &atomMap);
    bestAtLast = false;
    if (rms < bestRMS) {
      bestRMS = rms;
      bestMap = atomMap;
      bestAtLast = true;
    }
  }
  if (!bestAtLast) {
    alignMol(prbMol, refMol, prbCid, refCid, &bestMap);
  }
  return bestRMS;
}
}
}

/*


def GetConformerRMSMatrix(mol, atomIds=None, prealigned=False):
  """ Returns the RMS matrix of the conformers of a molecule.
  As a side-effect, the conformers will be aligned to the first
  conformer (i.e. the reference) and will left in the aligned state.

  Arguments:
    - mol:     the molecule
    - atomIds: (optional) list of atom ids to use a points for
               alingment - defaults to all atoms
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
      rmsvals.append(GetConformerRMS(mol, 0, i, atomIds=atomIds,
prealigned=prealigned))
  # loop over the conformations (except the reference one)
  cmat = []
  for i in range(1, mol.GetNumConformers()):
    cmat.append(rmsvals[i - 1])
    for j in range(1, i):
      cmat.append(GetConformerRMS(mol, i, j, atomIds=atomIds, prealigned=True))
  return cmat
*/
