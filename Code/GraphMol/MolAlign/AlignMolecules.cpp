// $Id$
//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AlignMolecules.h"
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace MolAlign {

// New API (uses AlignmentParameters)
double alignMol(ROMol &prbMol, const ROMol &refMol,
                const AlignmentParameters &alignParameter) {
  //  alignMolConformers(prbMol, alignParameter);
  RDGeom::Transform3D trans;
  double rmsd = getAlignmentTransform(prbMol, refMol, trans, alignParameter);

  // now transform the relevant conformation on prbMol
  Conformer &conf = prbMol.getConformer(alignParameter.prbConformerID);
  MolTransforms::transformConformer(conf, trans);
  return rmsd;
}

void alignMolConformers(ROMol &mol, const AlignmentParameters &alignParameter,
                        const std::vector<unsigned int> *atomIds,
                        const std::vector<unsigned int> *confIds,
                        std::vector<double> *RMSlist) {
  AlignmentParameters cAlignParameter = AlignmentParameters(alignParameter);
  cAlignParameter.enumerateAll = false;
  if (atomIds) {
    MatchVectType atomMap;
    BOOST_FOREACH (int atomId, *atomIds) {
      atomMap.push_back(std::pair<int, int>(atomId, atomId));
    }
    cAlignParameter.atomMap = &atomMap;
  }
  //  std::cout << '\n' << cAlignParameter;

  // Create list of confIds if missing
  std::vector<unsigned int> cConfIds;
  if (confIds) {
    cConfIds = *confIds;
  } else {
    for (unsigned int i = 0; i < mol.getNumConformers(); i++) {
      cConfIds.push_back(i);
    }
  }
  //  std::cout << "Conformations: ";
  //  printVector(std::cout, cConfIds);
  //  std::cout << '\n';
  cAlignParameter.refConformerID = cConfIds[0];
  BOOST_FOREACH (int confId, cConfIds) {
    //    std::cout << "Align conformation " << confId;
    if (confId == cAlignParameter.refConformerID) {
      //      std::cout << "skipped\n";
      continue;
    }
    cAlignParameter.prbConformerID = confId;
    double rmsd = alignMol(mol, mol, cAlignParameter);
    if (RMSlist) {
      RMSlist->push_back(rmsd);
    }
    //    std::cout << "aligned" << rmsd << "\n";
  }
}

void _getHeavyIndices(const ROMol &mol, std::set<int> &hAtoms) {
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    if ((*atomIt)->getAtomicNum() != 1) {
      hAtoms.insert((*atomIt)->getIdx());
    }
  }
}

void getAtomMappings(RWMol &refMol, RWMol &prbMol,
                     std::vector<MatchVectType> &atomMaps,
                     const AlignmentParameters &alignParameter) {
  atomMaps.clear();
  if (alignParameter.atomMap) {
    atomMaps.resize(1);
    if (alignParameter.ignoreHydrogens) {
      // Remove mappings between hydrogens from the user defined userAtomMap
      std::set<int> refIdx;
      std::set<int> prbIdx;
      _getHeavyIndices(refMol, refIdx);
      _getHeavyIndices(prbMol, prbIdx);
      MatchVectType atomMap;
      for (MatchVectType::const_iterator mapping =
               alignParameter.atomMap->begin();
           mapping != alignParameter.atomMap->end(); ++mapping) {
        if (prbIdx.count(mapping->first) == 1 and
            refIdx.count(mapping->second) == 1) {
          atomMap.push_back(
              std::pair<int, int>(mapping->first, mapping->second));
        }
      }
      atomMaps[0] = atomMap;
      MolOps::removeHs(refMol);
      MolOps::removeHs(prbMol);
    } else {
      // copy the userAtomMap as the first element in atomMaps
      atomMaps[0] = MatchVectType(*alignParameter.atomMap);
    }
  } else {
    if (alignParameter.ignoreHydrogens) {
      MolOps::removeHs(refMol);
      MolOps::removeHs(prbMol);
    }
    getSubstructureAtomMapping(refMol, prbMol, atomMaps, alignParameter);
  }
}

void getSubstructureAtomMapping(const ROMol &refMol, const ROMol &prbMol,
                                std::vector<MatchVectType> &atomMaps,
                                const AlignmentParameters &alignParameter) {
  const bool recursionPossible = true;
  const bool useChirality = false;
  const bool useQueryQueryMatches = true;
  atomMaps.clear();
  if (alignParameter.enumerateAll) {
    const bool uniquify = false;
    SubstructMatch(refMol, prbMol, atomMaps, uniquify, recursionPossible,
                   useChirality, useQueryQueryMatches);
  } else {
    MatchVectType match;
    SubstructMatch(refMol, prbMol, match, recursionPossible, useChirality,
                   useQueryQueryMatches);
    atomMaps.resize(1);
    atomMaps[0] = match;
  }
}

double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans, int prbCid, int refCid,
                             const MatchVectType &atomMap,
                             const AlignmentParameters &alignParameter) {
  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  const Conformer &prbCnf = prbMol.getConformer(prbCid);
  const Conformer &refCnf = refMol.getConformer(refCid);

  MatchVectType::const_iterator mi;
  for (mi = atomMap.begin(); mi != atomMap.end(); mi++) {
    prbPoints.push_back(&prbCnf.getAtomPos(mi->first));
    refPoints.push_back(&refCnf.getAtomPos(mi->second));
  }
  double ssr = RDNumeric::Alignments::AlignPoints(
      refPoints, prbPoints, trans, alignParameter.weights,
      alignParameter.reflect, alignParameter.maxIterations);
  ssr /= (prbPoints.size());
  return sqrt(ssr);
}

double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans,
                             const AlignmentParameters &alignParameter) {
  // Create a vector of possible atom mappings
  //  std::cout << alignParameter << '\n';
  std::vector<MatchVectType> atomMaps;
  RWMol cRefMol = RWMol(refMol);
  RWMol cPrbMol = RWMol(prbMol);
  getAtomMappings(cRefMol, cPrbMol, atomMaps, alignParameter);
  //  std::cout << "Number of atom maps : " << atomMaps.size() << '\n';
  //  std::cout << atomMaps;

  // Select on atom mapping
  double rmsd;
  if (alignParameter.findBestAtomMap) {
    MatchVectType bestMap;
    bool bestAtLast = false;
    rmsd = 1000000.0;
    BOOST_FOREACH (const MatchVectType &atomMap, atomMaps) {
      double rms = getAlignmentTransform(
          cPrbMol, cRefMol, trans, alignParameter.prbConformerID,
          alignParameter.refConformerID, atomMap, alignParameter);
      bestAtLast = false;
      if (rms < rmsd) {
        //        std::cout << rms << " " << rmsd << '\n';
        rmsd = rms;
        bestMap = atomMap;
        bestAtLast = true;
      }
    }
    if (!bestAtLast) {
      rmsd = getAlignmentTransform(
          cPrbMol, cRefMol, trans, alignParameter.prbConformerID,
          alignParameter.refConformerID, bestMap, alignParameter);
    }
  } else {
    rmsd = getAlignmentTransform(
        cPrbMol, cRefMol, trans, alignParameter.prbConformerID,
        alignParameter.refConformerID, atomMaps[0], alignParameter);
  }
  return rmsd;
}

double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans, int prbCid, int refCid,
                             const MatchVectType *atomMap,
                             const RDNumeric::DoubleVector *weights,
                             bool reflect, unsigned int maxIterations) {
  AlignmentParameters alignPara;
  alignPara.ignoreHydrogens = false;
  alignPara.enumerateAll = false;
  alignPara.atomMap = (const MatchVectType *)atomMap;
  alignPara.weights = weights;
  alignPara.reflect = reflect;
  alignPara.maxIterations = maxIterations;
  alignPara.prbConformerID = prbCid;
  alignPara.refConformerID = refCid;
  return getAlignmentTransform(prbMol, refMol, trans, alignPara);
}

double alignMol(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
                const MatchVectType *atomMap,
                const RDNumeric::DoubleVector *weights, bool reflect,
                unsigned int maxIterations) {
  AlignmentParameters alignPara;
  alignPara.ignoreHydrogens = false;
  alignPara.enumerateAll = false;
  alignPara.atomMap = (const MatchVectType *)atomMap;
  alignPara.weights = weights;
  alignPara.reflect = reflect;
  alignPara.maxIterations = maxIterations;
  alignPara.prbConformerID = prbCid;
  alignPara.refConformerID = refCid;
  return alignMol(prbMol, refMol, alignPara);
}

void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
                        const std::vector<unsigned int> *atomIds = 0) {
  unsigned int na = conf.getNumAtoms();
  pts.clear();
  if (atomIds == 0) {
    unsigned int ai;
    pts.reserve(na);
    for (ai = 0; ai < na; ++ai) {
      pts.push_back(&conf.getAtomPos(ai));
    }
  } else {
    pts.reserve(atomIds->size());
    std::vector<unsigned int>::const_iterator cai;
    for (cai = atomIds->begin(); cai != atomIds->end(); cai++) {
      pts.push_back(&conf.getAtomPos(*cai));
    }
  }
}

void alignMolConformers(ROMol &mol, const std::vector<unsigned int> *atomIds,
                        const std::vector<unsigned int> *confIds,
                        const RDNumeric::DoubleVector *weights, bool reflect,
                        unsigned int maxIters, std::vector<double> *RMSlist) {
  if (mol.getNumConformers() == 0) {
    // nothing to be done ;
    return;
  }
  AlignmentParameters alignPara;
  alignPara.ignoreHydrogens = false;
  alignPara.enumerateAll = false;
  alignPara.atomMap = NULL;
  alignPara.weights = weights;
  alignPara.reflect = reflect;
  alignPara.maxIterations = maxIters;
  alignMolConformers(mol, alignPara, atomIds, confIds, RMSlist);
}

void alignMolConformersOld(ROMol &mol, const std::vector<unsigned int> *atomIds,
                           const std::vector<unsigned int> *confIds,
                           const RDNumeric::DoubleVector *weights, bool reflect,
                           unsigned int maxIters,
                           std::vector<double> *RMSlist) {
  if (mol.getNumConformers() == 0) {
    // nothing to be done ;
    return;
  }

  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  int cid = -1;
  if ((confIds != 0) && (confIds->size() > 0)) {
    cid = confIds->front();
  }
  const Conformer &refCnf = mol.getConformer(cid);
  _fillAtomPositions(refPoints, refCnf, atomIds);

  // now loop through the remaining conformations and transform them
  RDGeom::Transform3D trans;
  double ssd;
  if (confIds == 0) {
    unsigned int i = 0;
    ROMol::ConformerIterator cnfi;
    // Conformer *conf;
    for (cnfi = mol.beginConformers(); cnfi != mol.endConformers(); cnfi++) {
      // conf = (*cnfi);
      i += 1;
      if (i == 1) {
        continue;
      }
      _fillAtomPositions(prbPoints, *(*cnfi), atomIds);
      ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
                                               weights, reflect, maxIters);
      if (RMSlist) {
        ssd /= (prbPoints.size());
        RMSlist->push_back(sqrt(ssd));
      }
      MolTransforms::transformConformer(*(*cnfi), trans);
    }
  } else {
    std::vector<unsigned int>::const_iterator cai;
    unsigned int i = 0;
    for (cai = confIds->begin(); cai != confIds->end(); cai++) {
      i += 1;
      if (i == 1) {
        continue;
      }
      Conformer &conf = mol.getConformer(*cai);
      _fillAtomPositions(prbPoints, conf, atomIds);
      ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
                                               weights, reflect, maxIters);
      if (RMSlist) {
        ssd /= (prbPoints.size());
        RMSlist->push_back(sqrt(ssd));
      }
      MolTransforms::transformConformer(conf, trans);
    }
  }
}

std::ostream &operator<<(std::ostream &os, const MatchVectType &map) {
  for (MatchVectType::const_iterator elem = map.begin(); elem != map.end();
       ++elem) {
    os << elem->first << '/' << elem->second << ' ';
  }
  os << '\n';
  return os;
}

std::ostream &operator<<(std::ostream &os,
                         const std::vector<MatchVectType> &maps) {
  BOOST_FOREACH (const MatchVectType &map, maps) { os << map; }
  return os;
}

void printVector(std::ostream &os, const std::vector<unsigned int> &idx) {
  BOOST_FOREACH (const unsigned int &id, idx) { os << id << ' '; }
}

std::ostream &operator<<(std::ostream &os, const AlignmentParameters &param) {
  os << "Alignment parameter\n";
  os << " findBestAtomMap " << param.findBestAtomMap << '\n';
  os << " ignoreHydrogens " << param.ignoreHydrogens << '\n';
  if (param.atomMap) {
    os << " atomMap " << *(param.atomMap);
  } else {
    os << " atomMap not defined\n";
  }
  return os;
}

}  // end namespace MolAlign
}
