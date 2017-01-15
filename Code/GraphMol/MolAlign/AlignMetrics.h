//
//  Copyright (C) 2001-2016 Peter Gedeck
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_ALIGNMETRICS_H_
#define _RD_ALIGNMETRICS_H_

#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <vector>

namespace RDKit {
namespace MolAlign {

//! Alignment metrics

//! Compute the optimal RMSD between two molecules
/*!
  Returns the optimal RMS for aligning two molecules, taking
  symmetry into account. As a side-effect, the probe molecule is
  left in the aligned state.

  \param refMol       the reference molecule
  \param prbMol       the molecule to be aligned to the reference
  \param refCid       (optional) reference conformation to use
  \param prbCid       (optional) probe conformation to use
  \param atomMaps     (optional) a list of atomMaps between the two molecules.
                      If not provided, these will be generated using a
                      substructure search.

  Note:
  This function will attempt to align all permutations of matching atom
  orders in both molecules, for some molecules it will lead to 'combinatorial
  explosion' especially if hydrogens are present.
  Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing the
  atom order.
  """
*/
double getBestRMS2(const ROMol &refMol, ROMol &prbMol, int probeConfId = -1,
                   int refConfId = -1, bool includeHydrogens = false,
                   std::vector<MatchVectType> *atomMaps = NULL);
}
}
#endif
