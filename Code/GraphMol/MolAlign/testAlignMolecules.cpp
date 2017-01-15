#include "AlignMolecules.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>

//#include <Geometry/point.h>
//#include <GraphMol/Atom.h>
//#include <GraphMol/Bond.h>
//#include <GraphMol/QueryBond.h>
//
//#include <GraphMol/FileParsers/MolSupplier.h>
//#include <GraphMol/FileParsers/MolWriters.h>
//#include <GraphMol/SmilesParse/SmilesWrite.h>
//#include <GraphMol/Descriptors/Crippen.h>
//#include <GraphMol/Conformer.h>
//#include <GraphMol/Substruct/SubstructMatch.h>
//#include <Numerics/Vector.h>
////#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
//
//#include <GraphMol/SmilesParse/SmilesParse.h>
//#include <GraphMol/MolTransforms/MolTransforms.h>

using namespace RDKit;

// Forward declaring a few utility functions
void compareConformers(const ROMol &mol, const Conformer &conf1,
                       const Conformer &conf2);
void get1oirAtomMap(MatchVectType &atomMap);

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

void test__getAlignmentTransform() {
  RWMol *refMol = SmilesToMol("CNC", 0, 1);
  RWMol *prbMol = SmilesToMol("CN", 0, 1);
  MolOps::addHs(*refMol, false, true);
  MolOps::addHs(*prbMol, false, true);

  INT_VECT refCids = DGeomHelpers::EmbedMultipleConfs(*refMol, 2);
  INT_VECT prbCids = DGeomHelpers::EmbedMultipleConfs(*prbMol, 2);

  MolAlign::AlignmentParameters alignPara;
  alignPara.enumerateAll = false;
  alignPara.ignoreHydrogens = true;

  double rmsd;
  RDGeom::Transform3D trans;
  std::vector<MatchVectType> matches;

  // We first look at aligning conformers of the same molecule
  MolAlign::getSubstructureAtomMapping(*refMol, *refMol, matches, alignPara);

  rmsd = MolAlign::_getAlignmentTransform(*refMol, *refMol, trans, 0, 0,
                                          matches[0], alignPara);
  TEST_ASSERT(rmsd == 0.0);
  rmsd = MolAlign::_getAlignmentTransform(*refMol, *refMol, trans, 0, 1,
                                          matches[0], alignPara);
  TEST_ASSERT(rmsd > 0.0);

  // Now align the two molecules (ignoring hydrogens)
  alignPara.ignoreHydrogens = true;
  MolAlign::getAtomMappings(*refMol, *prbMol, matches, alignPara);

  rmsd = MolAlign::_getAlignmentTransform(*prbMol, *refMol, trans, 0, 0,
                                          matches[0], alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.0, 0.01));
  rmsd = MolAlign::_getAlignmentTransform(*prbMol, *refMol, trans, 0, 1,
                                          matches[0], alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.0, 0.01));
}

void test_MolAlign() {
  std::string testDir = getenv("RDBASE");
  testDir += "/Code/GraphMol/MolAlign/test_data";
  ROMol *m1 = MolFileToMol(testDir + "/1oir.mol");
  ROMol *m2 = MolFileToMol(testDir + "/1oir_conf.mol");
  ROMol *m3 = MolFileToMol(testDir + "/1oir_trans.mol");

  MolAlign::AlignmentParameters alignPara;
  double rmsd = MolAlign::alignMol(*m2, *m1, alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // We can check the aligned conformer to 1oir_trans
  compareConformers(*m3, m2->getConformer(0), m3->getConformer(0));

  RDGeom::Transform3D trans;
  rmsd = MolAlign::getAlignmentTransform(*m1, *m2, trans);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // specify conformations
  alignPara.refConformerID = 0;
  alignPara.prbConformerID = 0;
  rmsd = MolAlign::alignMol(*m1, *m2, alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // Test the old API
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  delete m1;
  delete m2;
  delete m3;
}

void test_alignMolWithQuery() {
  // identical to test1MolAlign except we replace one atom with a QueryAtom
  // instead
  std::string testDir = getenv("RDBASE");
  testDir += "/Code/GraphMol/MolAlign/test_data";
  RWMol *m1 = MolFileToMol(testDir + "/1oir.mol");
  RWMol *m2 = MolFileToMol(testDir + "/1oir_conf.mol");
  RWMol *m3 = MolFileToMol(testDir + "/1oir_trans.mol");

  QueryAtom *a1 = new QueryAtom(6);
  QueryAtom *a2 = new QueryAtom(6);

  // we replace the same nitrogen instead with a null
  // query  28 and 19 are the "same" atoms
  m1->replaceAtom(28, a1);
  m2->replaceAtom(19, a2);
  m3->replaceAtom(0, new QueryAtom(5));

  MolAlign::AlignmentParameters alignPara;
  double rmsd = MolAlign::alignMol(*m2, *m1);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  compareConformers(*m3, m2->getConformer(0), m3->getConformer(0));

  RDGeom::Transform3D trans;
  rmsd = MolAlign::getAlignmentTransform(*m1, *m2, trans);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // specify conformations
  alignPara.refConformerID = 0;
  alignPara.prbConformerID = 0;
  rmsd = MolAlign::alignMol(*m1, *m2, alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // Test the old API
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  delete m1;
  delete m2;
  delete m3;
}

void test_alignMolAtomMap() {
  std::string testDir = getenv("RDBASE");
  testDir += "/Code/GraphMol/MolAlign/test_data";
  RWMol *m1 = MolFileToMol(testDir + "/1oir.mol");
  RWMol *m2 = MolFileToMol(testDir + "/1oir_conf.mol");

  MatchVectType atomMap;
  get1oirAtomMap(atomMap);

  MolAlign::AlignmentParameters alignPara;
  alignPara.atomMap = &atomMap;

  double rmsd = MolAlign::alignMol(*m2, *m1, alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.8525));

  // Test the old API
  rmsd = MolAlign::alignMol(*m2, *m1, 0, 0, &atomMap);
  TEST_ASSERT(RDKit::feq(rmsd, 0.8525));

  delete m1;
  delete m2;
}

void test_alignMolWeights() {
  std::string testDir = getenv("RDBASE");
  testDir += "/Code/GraphMol/MolAlign/test_data";
  RWMol *m1 = MolFileToMol(testDir + "/1oir.mol");
  RWMol *m2 = MolFileToMol(testDir + "/1oir_conf.mol");

  MatchVectType atomMap;
  get1oirAtomMap(atomMap);

  RDNumeric::DoubleVector wts(6);
  wts.setVal(0, 1.0);
  wts.setVal(1, 1.0);
  wts.setVal(2, 1.0);
  wts.setVal(3, 1.0);
  wts.setVal(4, 1.0);
  wts.setVal(5, 2.0);

  MolAlign::AlignmentParameters alignPara;
  alignPara.atomMap = &atomMap;
  alignPara.weights = &wts;

  double rmsd = MolAlign::alignMol(*m2, *m1, alignPara);
  TEST_ASSERT(RDKit::feq(rmsd, 0.9513));

  rmsd = MolAlign::alignMol(*m2, *m1, 0, 0, &atomMap, &wts);
  TEST_ASSERT(RDKit::feq(rmsd, 0.9513));

  delete m1;
  delete m2;
}

void test_alignMolIssue241() {
  std::string testDir = getenv("RDBASE");
  testDir += "/Code/GraphMol/MolAlign/test_data";
  RWMol *m1 = MolFileToMol(testDir + "/Issue241.mol");

  std::string res;
  MolPickler::pickleMol(*m1, res);
  ROMol *ref = new ROMol(res);
  DGeomHelpers::EmbedMolecule(*ref, 30, 239 * 10);
  ForceFields::ForceField *ff1 = UFF::constructForceField(*ref);
  ff1->initialize();
  ff1->minimize(200, 1e-8, 1e-6);

  std::string pkl2;
  MolPickler::pickleMol(*m1, pkl2);
  ROMol *probe = new ROMol(pkl2);
  DGeomHelpers::EmbedMolecule(*probe, 30, 239 * 10);
  ForceFields::ForceField *ff2 = UFF::constructForceField(*probe);
  ff2->initialize();
  ff2->minimize(200, 1e-8, 1e-6);

  double rmsd = MolAlign::alignMol(*ref, *probe);
  TEST_ASSERT(RDKit::feq(rmsd, 0.0));

  delete m1;
  delete ref;
  delete probe;
}

void compareConformers(const ROMol &mol, const Conformer &conf1,
                       const Conformer &conf2) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    RDGeom::Point3D pt1 = conf1.getAtomPos(i);
    RDGeom::Point3D pt2 = conf2.getAtomPos(i);
    TEST_ASSERT(RDKit::feq(pt1.x, pt2.x, 0.001));
    TEST_ASSERT(RDKit::feq(pt1.y, pt2.y, 0.001));
    TEST_ASSERT(RDKit::feq(pt1.z, pt2.z, 0.001));
  }
}

void get1oirAtomMap(MatchVectType &atomMap) {
  atomMap.clear();
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
}

void test_AlignMolecules() {
  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_getSubstructureAtomMapping \n\n";
  test_getSubstructureAtomMapping();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_getAtomMappings \n\n";
  test_getAtomMappings();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_getHeavyIndices \n\n";
  test_getHeavyIndices();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test__getAlignmentTransform \n\n";
  test__getAlignmentTransform();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_MolAlign \n\n";
  test_MolAlign();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_alignMolWithQuery \n\n";
  test_alignMolWithQuery();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_alignMolAtomMap \n\n";
  test_alignMolAtomMap();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_alignMolWeights \n\n";
  test_alignMolWeights();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test_alignMolIssue241 \n\n";
  test_alignMolIssue241();

  //  std::cout << "\t---------------------------------\n";
  //  std::cout << "\t test1AlignMetrics \n\n";
  //  test1AlignMetrics();
}
