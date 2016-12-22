# Namespace MolAlign

## Public functions/classes:

* class MolAlignException
* getAlignmentTransform: 
  * Compute the transformation required to align two molecules, probe (fixed) and reference (changed).
  * Returns RMSD and transformation matrix (trans)
* alignMol:
  * Optimally (minimum RMSD) align a molecule to another molecule
  * Returns RMSD
* alignMolConformers:
  * Align the conformations of a molecule using a common set of atoms.
  * Returns aligned conformations and list of RMSD values to reference conformation (optional)
* getConformerRMS:
  * Compute the RMSD between two conformations
  * Returns RMSD and aligns molecule
* getBestRMS:
  * Compute the optimal RMSD between two molecules

## Open3Dalign
* O3AFuncData:
  * Data structure that contains configuration options
* O3AConstraint:
  * Alignment constraints
* class MolHistogram:
* class LAP
* class SDM
* class O3A
* randomTransform
* o3aMMFFCostFunc
* o3aMMFFWeightFunc
* o3aMMFFScoringFunc
* o3aCrippenCostFunc
* o3aCrippenWeightFunc
* o3aCrippenScoringFunc
* getO3AForProbeConfs

