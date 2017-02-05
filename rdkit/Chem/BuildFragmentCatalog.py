#
#  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""  command line utility for working with FragmentCatalogs (CASE-type analysis)

**Usage**

  BuildFragmentCatalog [optional args] <filename>

 filename, the name of a delimited text file containing InData, is required
 for some modes of operation (see below)

**Command Line Arguments**

 - -n *maxNumMols*:  specify the maximum number of molecules to be processed

 - -b: build the catalog and OnBitLists
    *requires InData*

 - -s: score compounds
    *requires InData and a Catalog, can use OnBitLists*

 - -g: calculate info gains
    *requires Scores*

 - -d: show details about high-ranking fragments
    *requires a Catalog and Gains*

 - --catalog=*filename*: filename with the pickled catalog.
    If -b is provided, this file will be overwritten.

 - --onbits=*filename*: filename to hold the pickled OnBitLists.
   If -b is provided, this file will be overwritten

 - --scores=*filename*: filename to hold the text score data.
   If -s is provided, this file will be overwritten

 - --gains=*filename*: filename to hold the text gains data.
   If -g is provided, this file will be overwritten

 - --details=*filename*: filename to hold the text details data.
   If -d is provided, this file will be overwritten.

 - --minPath=2: specify the minimum length for a path

 - --maxPath=6: specify the maximum length for a path

 - --smiCol=1: specify which column in the input data file contains
     SMILES

 - --actCol=-1: specify which column in the input data file contains
     activities

 - --nActs=2: specify the number of possible activity values

 - --nBits=-1: specify the maximum number of bits to show details for

"""
from __future__ import print_function

import argparse
from contextlib import contextmanager, closing
from io import StringIO
import os
import sys
import time

import numpy

from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import InfoTheory
from rdkit.six import next
from rdkit.six.moves import cPickle as pickle

_MSGDEST = sys.stdout


def message(msg, dest=None):
  dest = dest or _MSGDEST
  dest.write(msg)


def BuildCatalog(suppl, maxPts=-1, groupFileName=None, minPath=2, maxPath=6, reportFreq=10):
  """ builds a fragment catalog from a set of molecules in a delimited text block

    **Arguments**

      - suppl: a mol supplier

      - maxPts: (optional) if provided, this will set an upper bound on the
        number of points to be considered

      - groupFileName: (optional) name of the file containing functional group
        information

      - minPath, maxPath: (optional) names of the minimum and maximum path lengths
        to be considered

      - reportFreq: (optional) how often to display status information

    **Returns**

      a FragmentCatalog

  """
  if groupFileName is None:
    groupFileName = os.path.join(RDConfig.RDDataDir, "FunctionalGroups.txt")

  fpParams = FragmentCatalog.FragCatParams(minPath, maxPath, groupFileName)
  catalog = FragmentCatalog.FragCatalog(fpParams)
  fgen = FragmentCatalog.FragCatGenerator()
  if maxPts > 0:
    nPts = maxPts
  else:
    if hasattr(suppl, '__len__'):
      nPts = len(suppl)
    else:
      nPts = -1
  progress = reportProgress(reportFreq, nPts)
  for i, mol in enumerate(suppl):
    if i == nPts:
      break
    progress('{0} paths'.format(catalog.GetFPLength()))
    fgen.AddFragsFromMol(mol, catalog)
  return catalog


def ScoreMolecules(suppl, catalog, maxPts=-1, actName='', acts=None, nActs=2, reportFreq=10):
  """ scores the compounds in a supplier using a catalog

    **Arguments**

      - suppl: a mol supplier

      - catalog: the FragmentCatalog

      - maxPts: (optional) the maximum number of molecules to be
        considered

      - actName: (optional) the name of the molecule's activity property.
        If this is not provided, the molecule's last property will be used.

      - acts: (optional) a sequence of activity values (integers).
        If not provided, the activities will be read from the molecules.

      - nActs: (optional) number of possible activity values

      - reportFreq: (optional) how often to display status information

    **Returns**

      a 2-tuple:

        1) the results table (a 3D array of ints nBits x 2 x nActs)

        2) a list containing the on bit lists for each molecule

  """
  nBits = catalog.GetFPLength()
  resTbl = numpy.zeros((nBits, 2, nActs), numpy.int)
  obls = []

  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]

  fpgen = FragmentCatalog.FragFPGenerator()
  suppl.reset()
  progress = reportProgress(reportFreq)
  for i, mol in enumerate(suppl, 1):
    progress()
    if mol:
      if not acts:
        act = int(mol.GetProp(actName))
      else:
        act = acts[i - 1]
      fp = fpgen.GetFPForMol(mol, catalog)
      bits = list(fp.GetOnBits())
      obls.append(bits)

      bits = [b - 1 for b in bits]
      resTbl[range(nBits), 0, act] += 1
      resTbl[bits, 0, act] -= 1
      resTbl[bits, 1, act] += 1
    else:
      obls.append([])
  return resTbl, obls


def ScoreFromLists(bitLists, suppl, catalog, maxPts=-1, actName='', acts=None, nActs=2,
                   reportFreq=10):
  """  similar to _ScoreMolecules()_, but uses pre-calculated bit lists
    for the molecules (this speeds things up a lot)


    **Arguments**

      - bitLists: sequence of on bit sequences for the input molecules

      - suppl: the input supplier (we read activities from here)

      - catalog: the FragmentCatalog

      - maxPts: (optional) the maximum number of molecules to be
        considered

      - actName: (optional) the name of the molecule's activity property.
        If this is not provided, the molecule's last property will be used.

      - nActs: (optional) number of possible activity values

      - reportFreq: (optional) how often to display status information

    **Returns**

       the results table (a 3D array of ints nBits x 2 x nActs)

  """
  nBits = catalog.GetFPLength()
  if maxPts > 0:
    nPts = maxPts
  else:
    nPts = len(bitLists)
  resTbl = numpy.zeros((nBits, 2, nActs), numpy.int)
  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]
  suppl.reset()

  progress = reportProgress(reportFreq, nPts, start=1)
  for i in range(1, nPts + 1):
    progress()
    mol = next(suppl)
    if not acts:
      act = int(mol.GetProp(actName))
    else:
      act = acts[i - 1]

    ids = list(set(ID - 1 for ID in bitLists[i - 1]))
    resTbl[range(nBits), 0, act] += 1
    resTbl[ids, 0, act] -= 1
    resTbl[ids, 1, act] += 1
  return resTbl


def CalcGains(suppl, catalog, topN=-1, actName='', acts=None, nActs=2, reportFreq=10, biasList=None,
              collectFps=0):
  """ calculates info gains by constructing fingerprints
    *DOC*

    Returns a 2-tuple:
       1) gains matrix
       2) list of fingerprints

  """
  nBits = catalog.GetFPLength()
  if topN < 0:
    topN = nBits
  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]

  if hasattr(suppl, '__len__'):
    nMols = len(suppl)
  else:
    nMols = -1
  fpgen = FragmentCatalog.FragFPGenerator()
  # ranker = InfoTheory.InfoBitRanker(nBits,nActs,InfoTheory.InfoType.ENTROPY)
  if biasList:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.BIASENTROPY)
    ranker.SetBiasList(biasList)
  else:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.ENTROPY)
  fps = []
  progress = reportProgress(reportFreq, nMols)
  for i, mol in enumerate(suppl, 0):
    if not acts:
      try:
        act = int(mol.GetProp(actName))
      except KeyError:
        message('ERROR: Molecule has no property: %s\n' % (actName))
        message('\tAvailable properties are: %s\n' % (str(mol.GetPropNames())))
        raise KeyError(actName)
    else:
      act = acts[i]
    progress()
    fp = fpgen.GetFPForMol(mol, catalog)
    ranker.AccumulateVotes(fp, act)
    if collectFps:
      fps.append(fp)
  gains = ranker.GetTopN(topN)
  return gains, fps


def CalcGainsFromFps(suppl, fps, topN=-1, actName='', acts=None, nActs=2, reportFreq=10,
                     biasList=None):
  """ calculates info gains from a set of fingerprints

    *DOC*

  """
  nBits = len(fps[0])
  if topN < 0:
    topN = nBits
  if not actName and not acts:
    actName = suppl[0].GetPropNames()[-1]

  if hasattr(suppl, '__len__'):
    nMols = len(suppl)
  else:
    nMols = -1
  if biasList:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.BIASENTROPY)
    ranker.SetBiasList(biasList)
  else:
    ranker = InfoTheory.InfoBitRanker(nBits, nActs, InfoTheory.InfoType.ENTROPY)
  progress = reportProgress(reportFreq, nMols)
  for i, mol in enumerate(suppl):
    if not acts:
      try:
        act = int(mol.GetProp(actName))
      except KeyError:
        message('ERROR: Molecule has no property: %s\n' % (actName))
        message('\tAvailable properties are: %s\n' % (str(mol.GetPropNames())))
        raise KeyError(actName)
    else:
      act = acts[i]
    progress()
    fp = fps[i]
    ranker.AccumulateVotes(fp, act)
  gains = ranker.GetTopN(topN)
  return gains


def reportProgress(reportFreq, nobjects=-1, start=0):
  """ Print progress messages """
  if nobjects > 0:
    fmt = 'Done {0} of {1}.\n'
    fmtMsg = 'Done {0} of {1}, {2}.\n'
  else:
    fmt = 'Done {0}.\n'
    fmtMsg = 'Done {0}, {2}.\n'
  count = start

  def counter(msg=None):
    nonlocal count
    if count and not count % reportFreq:
      if msg is None:
        message(fmt.format(count, nobjects))
      else:
        message(fmtMsg.format(count, nobjects, msg))
    count += 1

  return counter


def OutputGainsData(outF, gains, cat, nActs=2):
  actHeaders = ['Act-%d' % (x) for x in range(nActs)]
  if cat:
    outF.write('id,Description,Gain,%s\n' % (','.join(actHeaders)))
  else:
    outF.write('id,Gain,%s\n' % (','.join(actHeaders)))
  for entry in gains:
    id_ = int(entry[0])
    outL = [str(id_)]
    if cat:
      descr = cat.GetBitDescription(id_)
      outL.append(descr)
    outL.append('%.6f' % entry[1])
    outL += ['%d' % x for x in entry[2:]]
    outF.write(','.join(outL))
    outF.write('\n')


def ProcessGainsData(inF, delim=',', idCol=0, gainCol=1):
  """ reads a list of ids and info gains out of an input file

  """
  res = []
  header = inF.readline().strip().split(delim)
  gainCol = header.index('Gain')
  idCol = header.index('id')
  for line in inF:
    splitL = line.strip().split(delim)
    res.append((splitL[idCol], float(splitL[gainCol])))
  return res


def ShowDetails(catalog, gains, nToDo=-1, outF=sys.stdout, idCol=0, gainCol=1, outDelim=','):
  """
   gains should be a sequence of sequences.  The idCol entry of each
   sub-sequence should be a catalog ID.  _ProcessGainsData()_ provides
   suitable input.

  """
  if nToDo < 0:
    nToDo = len(gains)
  for i in range(nToDo):
    id_ = int(gains[i][idCol])
    gain = float(gains[i][gainCol])
    descr = catalog.GetBitDescription(id_)
    if descr:
      outF.write('%s\n' % (outDelim.join((str(id_), descr, str(gain)))))


def SupplierFromDetails(details):
  from rdkit.VLib.NodeLib.DbMolSupply import DbMolSupplyNode
  from rdkit.VLib.NodeLib.SmilesSupply import SmilesSupplyNode

  if details.dbName:
    conn = DbConnect(details.dbName, details.tableName)
    suppl = DbMolSupplyNode(conn.GetData())
  else:
    suppl = SmilesSupplyNode(details.inFileName, delim=details.delim, nameColumn=details.nameCol,
                             smilesColumn=details.smiCol, titleLine=details.hasTitle)
    if isinstance(details.actCol, int):
      suppl.reset()
      m = next(suppl)
      actName = m.GetPropNames()[details.actCol]
      details.actCol = actName
    if isinstance(details.nameCol, int):
      suppl.reset()
      m = next(suppl)
      nameName = m.GetPropNames()[details.nameCol]
      details.nameCol = nameName
      suppl.reset()
  if isinstance(details.actCol, int):
    suppl.reset()
    m = next(suppl)
    actName = m.GetPropNames()[details.actCol]
    details.actCol = actName
  if isinstance(details.nameCol, int):
    suppl.reset()
    m = next(suppl)
    nameName = m.GetPropNames()[details.nameCol]
    details.nameCol = nameName
    suppl.reset()
  return suppl

# def Usage():
#   print("This is BuildFragmentCatalog")
#   print('usage error')
#   # print(__doc__)
#   sys.exit(-1)


def initParser():
  """ Initialize the parser for the CLI """
  parser = argparse.ArgumentParser(
    description='command line utility for working with ' + 'FragmentCatalogs (CASE-type analysis)')
  parser.add_argument('-n', metavar='maxNumMols', type=int, default=-1, dest='numMols',
                      help='specify the maximum number of molecules to be processed')
  parser.add_argument('-d', default='', dest='dbName', help='Database name')
  group = parser.add_mutually_exclusive_group()
  group.add_argument('-c', default=',', dest='delim', action='store_const', const=',',
                     help='Read comma separated file')
  group.add_argument('-s', dest='delim', action='store_const', const=' ',
                     help='Read space separated file')
  group.add_argument('-t', dest='delim', action='store_const', const='\t',
                     help='Read tab separated file')
  parser.add_argument('--build', default=False, action='store_true', dest='doBuild',
                      help='build the catalog and OnBitLists (requires InData)')
  parser.add_argument('--sigs', default=False, action='store_true', dest='doSigs')
  parser.add_argument('--score', default=False, action='store_true', dest='doScore',
                      help='score compounds (requires InData and a Catalog, can use OnBitLists)')
  parser.add_argument('--gains', default=False, action='store_true', dest='doGains',
                      help='calculate info gains (requires Scores)')
  parser.add_argument(
    '--details', default=False, action='store_true', dest='doDetails',
    help='show details about high-ranking fragments ' + '(requires a Catalog and Gains)')
  parser.add_argument('--catalog', metavar='filename', dest='catalogName', default=None,
                      help='filename with the pickled catalog. ' +
                      'If --build is provided, this file will be overwritten.')
  parser.add_argument('--onbits', metavar='filename', dest='onBitsName', default=None,
                      help='filename to hold the pickled OnBitLists. ' +
                      'If -b is provided, this file will be overwritten.')
  parser.add_argument('--scoresFile', metavar='filename', dest='scoresName', default=None,
                      help='filename to hold the text score data. ' +
                      'If -s is provided, this file will be overwritten')
  parser.add_argument('--gainsFile', metavar='filename', dest='gainsName', default=None,
                      help='filename to hold the text gains data. ' +
                      'If -g is provided, this file will be overwritten')
  parser.add_argument('--detailsFile', metavar='filename', dest='detailsName', default=None,
                      help='filename to hold the text details data. ' +
                      'If -d is provided, this file will be overwritten.')
  parser.add_argument('--fpFile', metavar='filename', dest='fpName', default=None)
  parser.add_argument('--minPath', metavar='N', type=int, default=2,
                      help='specify the minimum length for a path')
  parser.add_argument('--maxPath', metavar='N', type=int, default=6,
                      help='specify the maximum length for a path')
  parser.add_argument(
    '--smiCol', metavar='N', default=1, type=intOrString,
    help='specify which column in the input data file contains SMILES (default %(default)s')
  parser.add_argument(
    '--actCol', metavar='N', default=-1, type=intOrString,
    help='specify which column in the input data file contains activities (default %(default)s')
  parser.add_argument('--nameCol', metavar='N', default=-1, type=intOrString)
  parser.add_argument('--nActs', metavar='N', default=2, type=int,
                      help='specify the number of possible activity values')
  parser.add_argument('--nBits', metavar='N', default=-1, type=int,
                      help='specify the maximum number of bits to show details for')
  parser.add_argument('--biasList', metavar='LIST', default=None, type=toTuple)
  parser.add_argument('--topN', metavar='N', default=-1, type=int)
  parser.add_argument('--noTitle', default=True, dest='hasTitle', action='store_false')
  parser.add_argument('input', default=None, help='File or table name')
  return parser


def intOrString(value):
  """ If possible convert the value to int, otherwise leave as string """
  try:
    return int(value)
  except ValueError:
    return value


def toTuple(value):
  """ Evaluate the value and assign to tuple """
  return tuple(eval(value))


def validateArgs(args, parser):
  """ Do some consistency checks """
  if args.dbName:
    args.tableName = args.input
  else:
    args.inFileName = args.input


@contextmanager
def timeIt():
  t1 = time.time()
  yield
  message("\tThat took %.2f seconds.\n" % (time.time() - t1))


def processDetails(details, parser):
  suppl = SupplierFromDetails(details)

  if details.doBuild or details.doScore or details.doGains:
    if not suppl:
      parser.error("We require inData to generate a catalog\n")

  cat = None
  obls = None
  if details.doBuild:
    message("Building catalog\n")
    with timeIt():
      cat = BuildCatalog(suppl, maxPts=details.numMols, minPath=details.minPath,
                         maxPath=details.maxPath)
    if details.catalogName:
      message("Dumping catalog data\n")
      with open(details.catalogName, 'wb+') as f:
        pickle.dump(cat, f)

  elif details.catalogName:
    message("Loading catalog\n")
    with open(details.catalogName, 'rb') as f:
      cat = pickle.load(f)
    if details.onBitsName:
      try:
        with open(details.onBitsName, 'rb') as f:
          obls = pickle.load(f)
      except Exception:
        obls = None
#       else:
#         if len(obls) < (inD.count('\n') - 1):
#           obls = None
  scores = None
  if details.doScore:
    if not cat:
      parser.error("We require a catalog to score molecules\n")
    message("Scoring compounds\n")
    if not obls or len(obls) < details.numMols:
      scores, obls = ScoreMolecules(suppl, cat, maxPts=details.numMols, actName=details.actCol,
                                    nActs=details.nActs)
      if details.scoresName:
        with open(details.scoresName, 'wb+') as f:
          pickle.dump(scores, f)
      if details.onBitsName:
        with open(details.onBitsName, 'wb+') as f:
          pickle.dump(obls, f)
    else:
      scores = ScoreFromLists(obls, suppl, cat, maxPts=details.numMols, actName=details.actCol,
                              nActs=details.nActs)
  elif details.scoresName:
    scores = pickle.load(open(details.scoresName, 'rb'))

  if details.fpName and os.path.exists(details.fpName) and not details.doSigs:
    message("Reading fingerprints from file.\n")
    with open(details.fpName, 'rb') as f:
      fps = pickle.load(f)
  else:
    fps = []
  gains = None
  if details.doGains:
    if not (cat or fps):
      parser.error("We require either a catalog or fingerprints to calculate gains\n")
    message("Calculating Gains\n")
    with timeIt():
      if details.fpName:
        collectFps = 1
      else:
        collectFps = 0
      if not fps:
        gains, fps = CalcGains(suppl, cat, topN=details.topN, actName=details.actCol,
                               nActs=details.nActs, biasList=details.biasList,
                               collectFps=collectFps)
        if details.fpName:
          message("Writing fingerprint file.\n")
          with open(details.fpName, 'wb+') as f:
            pickle.dump(fps, f, 1)
      else:
        gains = CalcGainsFromFps(suppl, fps, topN=details.topN, actName=details.actCol,
                                 nActs=details.nActs, biasList=details.biasList)
    if details.gainsName:
      with open(details.gainsName, 'w+') as f:
        OutputGainsData(f, gains, cat, nActs=details.nActs)
  else:
    if details.gainsName:
      with open(details.gainsName, 'r') as f:
        gains = ProcessGainsData(f)

  if details.doDetails:
    if not cat:
      parser.error("We require a catalog to get details\n")
    if gains is None:
      parser.error("We require gains data to get details\n")
    with closing(StringIO()) as io:
      io.write('id,SMILES,gain\n')
      ShowDetails(cat, gains, nToDo=details.nBits, outF=io)
      if details.detailsName:
        with open(details.detailsName, 'w+') as f:
          f.write(io.getvalue())
      else:
        sys.stderr.write(io.getvalue())

if __name__ == '__main__':
  parser = initParser()
  args = parser.parse_args()
  validateArgs(args, parser)
  processDetails(args, parser)
