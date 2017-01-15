import unittest
import doctest

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import MolAlign
import os

MolAlign.AlignmentParameters
MolAlign.AlignMolConformers
MolAlign.GetAlignmentTransform
MolAlign.GetBestRMS
MolAlign.GetConformerRMS

# # O3A alignment
# from rdkit.Chem.rdMolAlign import GetCrippenO3A, GetCrippenO3AForProbeConfs
# from rdkit.Chem.rdMolAlign import GetO3A, GetO3AForProbeConfs
# from rdkit.Chem.rdMolAlign import O3A, RandomTransform


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(MolAlign, optionflags=doctest.ELLIPSIS))
  return tests


class Test(unittest.TestCase):

  def testConformerRMS(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    _ = AllChem.EmbedMultipleConfs(m1, 2)

    m2 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    m2.AddConformer(m1.GetConformer(id=1))

    # test that the prealigned flag is working
    self.assertAlmostEqual(
      MolAlign.GetConformerRMS(m1, 0, 1, prealigned=True),
      MolAlign.pyGetConformerRMS(m1, 0, 1, prealigned=True))
    rms1 = MolAlign.GetConformerRMS(m1, 0, 1, prealigned=True)
    rms2 = MolAlign.GetConformerRMS(m1, 0, 1, prealigned=False)
    self.assertTrue(rms1 > rms2)
    self.assertAlmostEqual(
      MolAlign.GetConformerRMS(m1, 0, 1, prealigned=False),
      MolAlign.pyGetConformerRMS(m1, 0, 1, prealigned=False))

    # test that RMS is the same as calculated by AlignMol()
    self.assertAlmostEqual(rms2, MolAlign.GetBestRMS(m2, m1, 1, 0), 3)

    # the RMS with itself must be zero
    rms2 = MolAlign.GetConformerRMS(m1, 0, 0, prealigned=True)
    self.assertAlmostEqual(rms2, 0.0, 4)

    rms2 = MolAlign.GetBestRMS(m1, m1, 0, 0)
    self.assertAlmostEqual(rms2, 0.0, 4)

  def test_GetBestRMS(self):
    m1 = Chem.MolFromSmiles('CCC(C)(C)C')
    m1 = Chem.AddHs(m1)
    _ = AllChem.EmbedMultipleConfs(m1, 2)
    print(m1.GetNumAtoms())

    import time
    t0 = time.time()
    MolAlign.GetBestRMS(m1, m1, 1, 0)
    print(time.time() - t0)

    m1 = Chem.MolFromSmiles('CCC(C)(C)C')
    m1 = Chem.AddHs(m1)
    _ = AllChem.EmbedMultipleConfs(m1, 2)
    m2 = Chem.MolFromSmiles('CCC(C)(C)C')
    m2 = Chem.AddHs(m2)
    m2.AddConformer(m1.GetConformer(id=1))

    t0 = time.time()
    MolAlign.GetBestRMS(m1, m2, ignoreHydrogens=False)
    print(time.time() - t0)

  def testConformerRMSMatrix(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    _ = AllChem.EmbedMultipleConfs(m1, 3, randomSeed=0)

    # Check confIds
    rmsvals1 = []
    MolAlign.AlignMolConformers(m1, confIds=[2, 0, 1], RMSlist=rmsvals1)
    print(rmsvals1)
    rmsvals2 = []
    MolAlign.AlignMolConformers(m1, confIds=[0, 2, 1], RMSlist=rmsvals2)
    # aligning one conformer to the other can be reversed
    self.assertAlmostEqual(rmsvals1[0], rmsvals2[0], places=2)
    # but values are different if we align a conformer (1) to different conformers
    self.assertNotAlmostEqual(rmsvals1[1], rmsvals2[1], places=2)

    # Now we align everything to the first conformer before we calculate the full RMS matrix
    rmsvals = []
    MolAlign.AlignMolConformers(m1, RMSlist=rmsvals)
    #     print(rmsvals)
    #     print('old code ', [1.7259733364813168, 1.0050653394532882])

    rmat = MolAlign.GetConformerRMSMatrix(m1)
    self.assertEqual(len(rmat), 3)
    #     print(rmat)
    #     print('old code ', [1.7259733369398338, 1.0050653394532898, 2.1019488243436353])
    for rms1, rms2 in zip(rmsvals, rmat[:2]):
      self.assertAlmostEqual(rms1, rms2)

    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    _ = AllChem.EmbedMultipleConfs(m1, 3, randomSeed=0)

    rmatNoAlign = MolAlign.GetConformerRMSMatrix(m1, prealigned=True)
    rmatAligned = MolAlign.GetConformerRMSMatrix(m1, prealigned=False)
    for rms1, rms2 in zip(rmatNoAlign, rmatAligned):
      self.assertGreaterEqual(rms1, rms2)

    m1 = Chem.MolFromSmiles('NCCCC')
    _ = AllChem.EmbedMultipleConfs(m1, 3, randomSeed=0)
    # If we use only three atoms, the RMS values to the first conformer will be close to 0
    rmat = MolAlign.GetConformerRMSMatrix(m1, atomIds=[0, 1])
    self.assertAlmostEqual(rmat[0], rmat[1], places=5)

    return

    m2 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    m2.AddConformer(m1.GetConformer(id=0))

    # test that the RMS matrix has the correct size
    rmat = MolAlign.GetConformerRMSMatrix(m1)
    self.assertEqual(len(rmat), 3)
    print(rmat)
    print('old code ', [1.7259733369398338, 1.0050653394532898, 2.1019488243436353])

    # test the prealigned option
    rmat2 = MolAlign.GetConformerRMSMatrix(m1, prealigned=True)
    print(rmat2)
    self.assertAlmostEqual(rmat[0], rmat2[0])

    # test that the elements are in the right order
    self.assertAlmostEqual(rmat[0], MolAlign.GetBestRMS(m1, m2, 1, 0), 3)
    self.assertAlmostEqual(rmat[1], MolAlign.GetBestRMS(m1, m2, 2, 0), 3)

  def testGithub112(self):
    """ problems with AllChem.GetBestRMS() and molecules with Hs """
    testdata = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data')
    m0 = Chem.MolFromMolFile(os.path.join(testdata, 'github112_tgt.mol'), removeHs=False)
    m1 = Chem.MolFromMolFile(os.path.join(testdata, 'github112_qry.mol'), removeHs=False)
    alignPara = MolAlign.AlignmentParameters()
    alignPara.ignoreHydrogens = False
    rms = MolAlign.GetBestRMS(m0, m1, ignoreHydrogens=False)
    self.assertAlmostEqual(rms, 0.456, 3)
    rms = MolAlign.GetBestRMS(m0, m1)
    self.assertAlmostEqual(rms, 0.280, 3)


if __name__ == "__main__":  # pragma: nocover
  unittest.main()
