#include "AlignMolecules.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

void test_getAtomMappings() {
  std::vector<MatchVectType> matches;
  bool enumerateAll, ignoreH;
  RWMol *refMol = SmilesToMol("CNC", 0, 1);
  RWMol *prbMol = SmilesToMol("CN", 0, 1);
  MolOps::addHs(*refMol, false, true);
  MolOps::addHs(*prbMol, false, true);

  MolAlign::AlignmentParameters alignPara;

  // Auto-generation of all mappings based on substructure mappings
  alignPara.enumerateAll = true;
  enumerateAll = true;

  // We consider all atoms - for different molecules, this will create zero
  // matches
  alignPara.ignoreHydrogens = false;
  ignoreH = false;
  MolAlign::getAtomMappings(*refMol, *prbMol, matches, alignPara);
  TEST_ASSERT(refMol->getNumAtoms() == 10);
  TEST_ASSERT(matches.size() == 0);

  MolAlign::getAtomMappings(*refMol, *refMol, matches, alignPara);
  TEST_ASSERT(refMol->getNumAtoms() == 10);
  TEST_ASSERT(matches.size() == 72);

  // We consider heavy atoms only
  alignPara.ignoreHydrogens = true;
  MolAlign::getAtomMappings(*refMol, *prbMol, matches, alignPara);
  TEST_ASSERT(refMol->getNumAtoms() == 3);
  TEST_ASSERT(prbMol->getNumAtoms() == 2);
  TEST_ASSERT(matches.size() == 2);

  MolAlign::getAtomMappings(*refMol, *refMol, matches, alignPara);
  TEST_ASSERT(matches.size() == 2);
  delete refMol;
  delete prbMol;

  // Check with user defined atom mapping
  refMol = SmilesToMol("CNC", 0, 1);
  prbMol = SmilesToMol("CN", 0, 1);
  MolOps::addHs(*refMol, false, true);
  MolOps::addHs(*prbMol, false, true);

  MatchVectType userAtomMap;
  userAtomMap.resize(4);
  userAtomMap[0] = std::pair<int, int>(0, 1);  // heavy-heavy
  userAtomMap[1] = std::pair<int, int>(0, 4);  // heavy-H
  userAtomMap[2] = std::pair<int, int>(5, 0);  // H-heavy
  userAtomMap[3] = std::pair<int, int>(4, 6);  // H-H

  alignPara.atomMap = &userAtomMap;
  alignPara.ignoreHydrogens = false;
  MolAlign::getAtomMappings(*refMol, *prbMol, matches, alignPara);
  TEST_ASSERT(refMol->getNumAtoms() == 10);
  TEST_ASSERT(prbMol->getNumAtoms() == 7);
  TEST_ASSERT(matches.size() == 1);
  TEST_ASSERT(matches[0] == userAtomMap);
  TEST_ASSERT(matches[0].size() == 4);

  alignPara.ignoreHydrogens = true;
  MolAlign::getAtomMappings(*refMol, *prbMol, matches, alignPara);
  TEST_ASSERT(refMol->getNumAtoms() == 3);
  TEST_ASSERT(prbMol->getNumAtoms() == 2);
  // We only keep matches that involve heavy atoms
  TEST_ASSERT(matches.size() == 1);
  TEST_ASSERT(matches[0].size() == 1)
  TEST_ASSERT(matches[0][0] == userAtomMap[0]);
}

void test_getSubstructureAtomMapping() {
  ROMol *refMol = SmilesToMol("CNC", 0, 1);
  ROMol *prbMol = SmilesToMol("CN", 0, 1);
  std::vector<MatchVectType> matches;
  MolAlign::AlignmentParameters alignPara;
  alignPara.enumerateAll = false;
  MolAlign::getSubstructureAtomMapping(*refMol, *prbMol, matches, alignPara);
  TEST_ASSERT(matches.size() == 1)
  alignPara.enumerateAll = true;
  MolAlign::getSubstructureAtomMapping(*refMol, *prbMol, matches, alignPara);
  TEST_ASSERT(matches.size() == 2)

  // Now add the hydrogens and see the mapping explode
  ROMol *refMolH = MolOps::addHs(*refMol, false, true);
  ROMol *prbMolH = MolOps::addHs(*prbMol, false, true);
  alignPara.enumerateAll = false;
  MolAlign::getSubstructureAtomMapping(*refMolH, *prbMolH, matches, alignPara);
  TEST_ASSERT(matches.size() == 1)
  alignPara.enumerateAll = true;
  MolAlign::getSubstructureAtomMapping(*refMolH, *prbMolH, matches, alignPara);
  TEST_ASSERT(matches.size() == 0)
  alignPara.enumerateAll = false;
  MolAlign::getSubstructureAtomMapping(*refMolH, *refMolH, matches, alignPara);
  TEST_ASSERT(matches.size() == 1)
  alignPara.enumerateAll = true;
  MolAlign::getSubstructureAtomMapping(*refMolH, *refMolH, matches, alignPara);
  TEST_ASSERT(matches.size() == 72)

  delete refMol;
  delete prbMol;
  delete refMolH;
  delete prbMolH;
}

void test_getHeavyIndices() {
  const ROMol *mol = SmilesToMol("CN", 0, 1);
  const ROMol *molH = MolOps::addHs(*mol, false, true);

  std::set<int> idx;
  MolAlign::_getHeavyIndices(*mol, idx);
  TEST_ASSERT(idx.size() == 2);

  for (ROMol::ConstAtomIterator atomIt = (*mol).beginAtoms();
       atomIt != (*mol).endAtoms(); ++atomIt) {
    if (idx.count((*atomIt)->getIdx())) {
      TEST_ASSERT((*atomIt)->getAtomicNum() != 1);
    } else {
      TEST_ASSERT((*atomIt)->getAtomicNum() == 1);
    }
  }
}
