//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_ALIGNMOLECULES_H_
#define _RD_ALIGNMOLECULES_H_

#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <vector>

namespace RDKit {
class Conformer;
class ROMol;
class RWMol;

typedef std::vector<std::pair<int, int> > MatchVectType;

// Utility function for debug output
std::ostream &operator<<(std::ostream &os, const MatchVectType &map);
std::ostream &operator<<(std::ostream &os,
                         const std::vector<MatchVectType> &maps);
namespace MolAlign {

void printVector(std::ostream &os, const std::vector<unsigned int> &idx);

class MolAlignException : public std::exception {
 public:
  //! construct with an error message
  MolAlignException(const char *msg) : _msg(msg){};
  //! construct with an error message
  MolAlignException(const std::string msg) : _msg(msg){};
  //! get the error message
  const char *message() const { return _msg.c_str(); };
  ~MolAlignException() throw(){};

 private:
  std::string _msg;
};

//! Parameter object for controlling molecule alignments
/*!
  findBestAtomMap Try all atom maps and use the mapping with the smallest RMSD
                  for the alignment
  ignoreHydrogens Include hydrogens in substructure atom mapping
  enumerateAll    Enumerate all substructure mappings
  atomMap         A vector of pairs of atom IDs (probe AtomId, reference AtomId)
                  used to compute the alignments. If this mapping is
                  not specified an attempt is made to generate one by
                  substructure matching

  weights         Optionally specify weights for each of the atom pairs
  reflect         If true reflect the conformation of the probe molecule
  maxIterations   Maximum number of iteration used in minimizing the RMSD

  refConformerID   ID of reference conformer (default -1, to consider the first)
  prbConformerID   ID of probe conformer (default -1, to consider the first)
  refAllConformers Use all conformers of reference to find the best alignment
  prbAllConformers Use all conformers of probe to find the best alignment
*/
struct AlignmentParameters {
  bool findBestAtomMap;
  bool ignoreHydrogens;
  bool enumerateAll;
  const MatchVectType *atomMap;
  const RDNumeric::DoubleVector *weights;
  bool reflect;
  unsigned int maxIterations;
  int refConformerID;
  int prbConformerID;
  bool refAllConformers;
  bool prbAllConformers;
  AlignmentParameters()
      : findBestAtomMap(false),
        ignoreHydrogens(false),
        enumerateAll(false),
        atomMap(NULL),
        weights(NULL),
        reflect(false),
        maxIterations(50),
        refConformerID(-1),
        prbConformerID(-1),
        refAllConformers(false),
        prbAllConformers(false){};
};

std::ostream &operator<<(std::ostream &os, const AlignmentParameters &param);

/*
 * New API (uses AlignmentParameters)
 */

//! Alignment functions

//! Optimally (minimum RMSD) align a molecule to another molecule
/*!
  The 3D transformation required to align the specified conformation (default
  first) in the probe molecule to a specified conformation in the reference
  molecule is computed so that the root mean squared distance between a
  specified set of atoms is minimized. This transforms is them applied to the
  specified conformation (default first) in the probe molecule.

  \param prbMol          molecule that is to be aligned
  \param refMol          molecule used as the reference for the alignment
  \param alignParameter  control parameter for alignment process

  <b>Returns</b>
  RMSD value
*/
double alignMol(ROMol &prbMol, const ROMol &refMol,
                const AlignmentParameters &alignParameter);

//! Align the conformations of a molecule using a common set of atoms.
/*!
  If the molecules contains queries, then the queries must also match exactly.

  \param mol                 the molecule of interest.
  \param alignmentParameter  control parameter for alignment process
  \param atomIds             vector of atoms to be used to generate the
                             alignment - defaults to all
  \param confIds             vector of conformations to align - defaults to all
  \param RMSlist             if nonzero, this will be used to return the RMS
                             values between the reference conformation and the
                             other aligned conformations
*/
void alignMolConformers(ROMol &mol, const AlignmentParameters &alignParameter,
                        const std::vector<unsigned int> *atomIds = NULL,
                        const std::vector<unsigned int> *confIds = NULL,
                        std::vector<double> *RMSlist = 0);

//! Compute the RMSD between two conformations
/*!
  Returns the RMS between two conformations.
  By default, the conformers will be aligned to the first conformer
  of the molecule (i.e. the reference) before RMS calculation and,
  as a side-effect, will be left in the aligned state.

  \param mol        molecule with several conformations
  \param confId1    the id of the first conformer
  \param confId2    the id of the second conformer
  \param atomIds    list of atom ids to use as points for alignment (defaults
                    to all atoms)
  \param prealigned by default the conformers are assumed to be unaligned and
                    will therefore be aligned to the first conformer

  <b>Returns</b>
  RMSD value
*/
double getConformerRMS(ROMol &mol, unsigned int confId1, unsigned int confId2,
                       const std::vector<unsigned int> *atomIds = NULL,
                       bool prealigned = false);

//! Compute the transformation required to align a molecule
/*!
  The 3D transformation required to align the specied conformation in the probe
  molecule (alignParameter.prbConformerID) to a specified conformation in the
  reference molecule (alignParameter.refConformerID) is computed so that the
  root mean squared distance between a specified set of atoms is minimized.

  If alignParameter.findBestAtomMap is true, all atom maps (user provided or
  automatically generated) are tried and the one with the lowest RMSD used for
  the 3D transformation. Note that in this case it is best to ignore hydrogens
  to reduce the number of possible maps.

  \param prbMol          molecule that is to be aligned
  \param refMol          molecule used as the reference for the alignment
  \param trans           storage for the computed transform
  \param alignParameter  control parameter for alignment process

  <b>Returns</b>
  RMSD value
*/
double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans,
                             const AlignmentParameters &alignParameter);

/*
 * Old API
 */
//! Compute the transformation required to align a molecule
/*!
  The 3D transformation required to align the specied conformation in the probe
  molecule
  to a specified conformation in the reference molecule is computed so that the
  root mean
  squared distance between a specified set of atoms is minimized

  \param prbMol    molecule that is to be aligned
  \param refMol    molecule used as the reference for the alignment
  \param trans     storage for the computed transform
  \param prbCid    ID of the conformation in the probe to be used
                   for the alignment (defaults to first conformation)
  \param refCid    ID of the conformation in the ref molecule to which
                   the alignment is computed (defaults to first conformation)
  \param atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                   used to compute the alignments. If this mapping is
                   not specified an attempt is made to generate on by
                   substructure matching
  \param weights   Optionally specify weights for each of the atom pairs
  \param reflect   if true reflect the conformation of the probe molecule
  \param maxIters  maximum number of iteration used in mimizing the RMSD

  <b>Returns</b>
  RMSD value
*/
double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans, int prbCid = -1,
                             int refCid = -1, const MatchVectType *atomMap = 0,
                             const RDNumeric::DoubleVector *weights = 0,
                             bool reflect = false, unsigned int maxIters = 50);

//! Optimally (minimum RMSD) align a molecule to another molecule
/*!
  The 3D transformation required to align the specied conformation in the probe
  molecule
  to a specified conformation in the reference molecule is computed so that the
  root mean
  squared distance between a specified set of atoms is minimized. This
  transforms is them
  applied to the specified conformation in the probe molecule

  \param prbMol    molecule that is to be aligned
  \param refMol    molecule used as the reference for the alignment
  \param prbCid    ID of the conformation in the probe to be used
                   for the alignment (defaults to first conformation)
  \param refCid    ID of the conformation in the ref molecule to which
                   the alignment is computed (defaults to first conformation)
  \param atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                   used to compute the alignments. If this mapping is
                   not specified an attempt is made to generate on by
                   substructure matching
  \param weights   Optionally specify weights for each of the atom pairs
  \param reflect   if true reflect the conformation of the probe molecule
  \param maxIters  maximum number of iteration used in mimizing the RMSD

  <b>Returns</b>
  RMSD value
*/
double alignMol(ROMol &prbMol, const ROMol &refMol, int prbCid = -1,
                int refCid = -1, const MatchVectType *atomMap = 0,
                const RDNumeric::DoubleVector *weights = 0,
                bool reflect = false, unsigned int maxIters = 50);

//! Align the conformations of a molecule using a common set of atoms. If
// the molecules contains queries, then the queries must also match exactly.
/*!
  \param mol       The molecule of interest.
  \param atomIds   vector of atoms to be used to generate the alignment.
                   All atoms will be used is not specified
  \param confIds   vector of conformations to align - defaults to all
  \param weights   vector of weights to applied to particular atom pairs
                   defaults to all weights = 1
  \param reflect   toggles reflecting (about the origin) the alignment
  \param maxIters  the maximum number of iterations to attempt
  \param RMSlist   if nonzero, this will be used to return the RMS values
                   between the reference conformation and the other aligned
                   conformations
*/
void alignMolConformers(ROMol &mol,
                        const std::vector<unsigned int> *atomIds = 0,
                        const std::vector<unsigned int> *confIds = 0,
                        const RDNumeric::DoubleVector *weights = 0,
                        bool reflect = false, unsigned int maxIters = 50,
                        std::vector<double> *RMSlist = 0);

// private functions

//! Create one or more mappings between the two molecules
//! (alignParameter.enumerateAll). prbMol can be a substructure of refMol.
/*!
   \param refMol
   \param prbMol
   \param atomMaps
   \param alignParameter
 */
void getSubstructureAtomMapping(const ROMol &refMol, const ROMol &prbMol,
                                std::vector<MatchVectType> &atomMaps,
                                const AlignmentParameters &alignParameter);

//! Create one or more mappings between the two molecules.
/* The mapping is controlled by alignParameter. The alignParameter.atomMap
 allows the user to provide a mapping between the two molecules. Otherwise
 substructure mappings of prbMol in refMol are determined (use
 alignParameter.enumerateAll to get all possible mappings). The parameter
 alignParameter.ignoreHydrogens can be used to reduce the mappings to heavy
 atoms only.

   \param refMol
   \param prbMol
   \param atomMaps
   \param alignParameter
 */
void getAtomMappings(RWMol &refMol, RWMol &prbMol,
                     std::vector<MatchVectType> &atomMaps,
                     const AlignmentParameters &alignParameter);

// void alignMolConformers(
//    ROMol &mol, const std::vector<unsigned int> *atomIds = NULL,
//    const std::vector<unsigned int> *confIds = NULL,
//    const AlignmentParameters &alignParameter = AlignmentParameters(),
//    std::vector<double> *RMSlist = NULL);

// Private methods
void _getHeavyIndices(const ROMol &mol, std::set<int> &heavyAtoms);
double _getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                              RDGeom::Transform3D &trans, int prbCid,
                              int refCid, const MatchVectType &atomMap,
                              const AlignmentParameters &alignParameter);

}  // end namespace MolAlign
}
#endif
