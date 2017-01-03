# Alignment code

## Requirements
* Atoms considered:
  * All atoms
  * Heavy atoms only
* Atom matching:
  * User defined (std::vector<MatchVectType> or MatchVectType)
  * Single substructure match (MatchVectType)
  * Multiple substruture matches ( std::vector<MatchVectType>)
* Alignments (ref <-> prb):
  * Conformer <-> Conformer
  * Conformer <-> All/range conformers

## Relevant functions in other namespaces
* AlignPoints (RDNumeric::Alignments)
  * Compute an optimal alignment (minimum sum of squared distance) between two sets of points in 3D
  * Returns SSR (RMSE = sqrt(SSR/NUMPOINTS) and transformation matrix (trans)
* MolTransforms::transformConformer
  * Transform the conformation using the specified transformation
void transformConformer(RDKit::Conformer &conf,


## Namespace MolAlign
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

# Open3Dalign
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

