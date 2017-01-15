# $Id$
#
"""unit testing code for 3D stuff

"""
import os
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import TorsionFingerprints


class TestCase(unittest.TestCase):

  def test_MatchFunctions(self):
    mol = Chem.MolFromSmiles('CCCCCC')

    inv = [1, 1, 1, 2, 2, 3]
    # All different
    atoms = [mol.GetAtomWithIdx(idx) for idx in (0, 3, 5)]
    self.assertEqual(TorsionFingerprints._doMatch(inv, atoms), False)
    self.assertEqual(TorsionFingerprints._doNotMatch(inv, atoms), True)
    self.assertEqual(TorsionFingerprints._doMatchExcept1(inv, atoms), None)

    # All the same
    atoms = [mol.GetAtomWithIdx(idx) for idx in (0, 1, 2)]
    self.assertEqual(TorsionFingerprints._doMatch(inv, atoms), True)
    self.assertEqual(TorsionFingerprints._doNotMatch(inv, atoms), False)
    self.assertEqual(TorsionFingerprints._doMatchExcept1(inv, atoms), None)

    # One is different
    atoms = [mol.GetAtomWithIdx(idx) for idx in (0, 1, 3)]
    self.assertEqual(TorsionFingerprints._doMatch(inv, atoms), False)
    self.assertEqual(TorsionFingerprints._doNotMatch(inv, atoms), False)
    self.assertEqual(TorsionFingerprints._doMatchExcept1(inv, atoms), atoms[2])
    atoms = [mol.GetAtomWithIdx(idx) for idx in (0, 3, 1)]
    self.assertEqual(TorsionFingerprints._doMatchExcept1(inv, atoms), atoms[1])
    atoms = [mol.GetAtomWithIdx(idx) for idx in (3, 0, 1)]
    self.assertEqual(TorsionFingerprints._doMatchExcept1(inv, atoms), atoms[0])
    self.assertRaises(ValueError, TorsionFingerprints._doMatchExcept1, inv, [1, 2, 3, 4])

  def testTorsionFingerprints(self):
    # we use the xray structure from the paper (JCIM, 52, 1499, 2012): 1DWD
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1DWD_ligand.pdb')
    ref = Chem.MolFromSmiles(
      'NC(=[NH2+])c1ccc(C[C@@H](NC(=O)CNS(=O)(=O)c2ccc3ccccc3c2)C(=O)N2CCCCC2)cc1')
    mol = Chem.MolFromPDBFile(refFile)
    mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)

    # the torsion lists
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    self.assertEqual(len(tors_list), 11)
    self.assertEqual(len(tors_list_rings), 4)
    self.assertAlmostEqual(tors_list[-1][1], 180.0, 4)
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol, maxDev='spec')
    self.assertAlmostEqual(tors_list[-1][1], 90.0, 4)
    self.assertRaises(ValueError, TorsionFingerprints.CalculateTorsionLists, mol, maxDev='test')
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol, symmRadius=0)
    self.assertEqual(len(tors_list[0][0]), 2)

    # the weights
    weights = TorsionFingerprints.CalculateTorsionWeights(mol)
    self.assertAlmostEqual(weights[4], 1.0)
    self.assertEqual(len(weights), len(tors_list + tors_list_rings))
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, 15, 14)
    self.assertAlmostEqual(weights[3], 1.0)
    self.assertRaises(ValueError, TorsionFingerprints.CalculateTorsionWeights, mol, 15, 3)

    # the torsion angles
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    torsions = TorsionFingerprints.CalculateTorsionAngles(mol, tors_list, tors_list_rings)
    self.assertEqual(len(weights), len(torsions))
    self.assertAlmostEqual(torsions[2][0][0], 232.5346, 4)

    # the torsion fingerprint deviation
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions)
    self.assertAlmostEqual(tfd, 0.0)
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1PPC_ligand.pdb')
    mol2 = Chem.MolFromPDBFile(refFile)
    mol2 = AllChem.AssignBondOrdersFromTemplate(ref, mol2)
    torsions2 = TorsionFingerprints.CalculateTorsionAngles(mol2, tors_list, tors_list_rings)
    weights = TorsionFingerprints.CalculateTorsionWeights(mol)
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions2, weights=weights)
    self.assertAlmostEqual(tfd, 0.0691, 4)
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions2)
    self.assertAlmostEqual(tfd, 0.1115, 4)

    # the wrapper functions
    tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol, mol2)
    self.assertAlmostEqual(tfd, 0.0691, 4)

    mol.AddConformer(mol2.GetConformer(), assignId=True)
    mol.AddConformer(mol2.GetConformer(), assignId=True)
    tfd = TorsionFingerprints.GetTFDBetweenConformers(mol, confIds1=[0], confIds2=[1, 2])
    self.assertEqual(len(tfd), 2)
    self.assertAlmostEqual(tfd[0], 0.0691, 4)

    tfdmat = TorsionFingerprints.GetTFDMatrix(mol)
    self.assertEqual(len(tfdmat), 3)

  def testTorsionFingerprintsAtomReordering(self):
    # we use the xray structure from the paper (JCIM, 52, 1499, 2012): 1DWD
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1DWD_ligand.pdb')
    ref = Chem.MolFromSmiles(
      'NC(=[NH2+])c1ccc(C[C@@H](NC(=O)CNS(=O)(=O)c2ccc3ccccc3c2)C(=O)N2CCCCC2)cc1')
    mol1 = Chem.MolFromPDBFile(refFile)
    mol1 = AllChem.AssignBondOrdersFromTemplate(ref, mol1)

    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1DWD_ligand_reordered.pdb')
    mol2 = Chem.MolFromPDBFile(refFile)
    mol2 = AllChem.AssignBondOrdersFromTemplate(ref, mol2)

    tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol1, mol2)
    self.assertEqual(tfd, 0.0)

  def testTorsionFingerprintsColinearBonds(self):
    # test that single bonds adjacent to triple bonds are ignored
    mol = Chem.MolFromSmiles('CCC#CCC')
    tors_list, _ = TorsionFingerprints.CalculateTorsionLists(mol, ignoreColinearBonds=True)
    self.assertEqual(len(tors_list), 0)
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, ignoreColinearBonds=True)
    self.assertEqual(len(weights), 0)

    # test that they are not ignored, but alternative atoms searched for
    tors_list, _ = TorsionFingerprints.CalculateTorsionLists(mol, ignoreColinearBonds=False)
    self.assertEqual(len(tors_list), 1)
    self.assertEqual(tors_list[0][0][0], (0, 1, 4, 5))
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, ignoreColinearBonds=False)
    self.assertEqual(len(weights), 1)

    # test that single bonds adjacent to terminal triple bonds are always ignored
    mol = Chem.MolFromSmiles('C#CCC')
    tors_list, _ = TorsionFingerprints.CalculateTorsionLists(mol, ignoreColinearBonds=True)
    self.assertEqual(len(tors_list), 0)
    tors_list, _ = TorsionFingerprints.CalculateTorsionLists(mol, ignoreColinearBonds=False)
    self.assertEqual(len(tors_list), 0)

  def assertBondStereoRoundTrips(self, fname):
    path = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', fname)
    mol = Chem.MolFromMolFile(path)
    refSmiles = mol.GetProp("_Name")
    self.assertTrue(len(refSmiles) > 0)
    self.assertEqual(Chem.MolToSmiles(mol, isomericSmiles=True), refSmiles)

    # now test Chem.DetectBondStereoChemistry more directly by constructing the
    # molecule from scratch
    oldconf = mol.GetConformer(0)
    newconf = Chem.Conformer(mol.GetNumAtoms())
    newmol = Chem.RWMol()

    for atm in mol.GetAtoms():
      ratm = Chem.Atom(atm.GetAtomicNum())
      ratm.SetFormalCharge(atm.GetFormalCharge())
      newmol.AddAtom(ratm)

      atomidx = atm.GetIdx()
      pos = oldconf.GetAtomPosition(atomidx)
      newconf.SetAtomPosition(atomidx, pos)

    for bnd in mol.GetBonds():
      newmol.AddBond(bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx(), Chem.BondType(bnd.GetBondType()))
    newmol.AddConformer(newconf)

    Chem.SanitizeMol(newmol)
    Chem.DetectBondStereoChemistry(newmol, newconf)

    # these aren't necessary for this specific test case, but are for
    # a more general conversion routine, so would like to see them
    # tested eventually
    # Chem.AssignAtomChiralTagsFromStructure(newmol)
    # Chem.AssignStereochemistry(newmol)

    self.assertEqual(Chem.MolToSmiles(newmol, isomericSmiles=True), refSmiles)

  def testDetectBondStereoChemistry(self):
    self.assertBondStereoRoundTrips('cis.sdf')
    self.assertBondStereoRoundTrips('trans.sdf')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
