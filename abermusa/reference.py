''' Encapsulate a single alignment reference protein structure. '''

# Date:   Jun 25 2015
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from labblouin import IO
from cPickle import load

class reference:

  def __init__(self,protein_file_path):
    self.protein_path = protein_file_path
    self.folder_name = IO.makeFolderFromFileName(protein_file_path)

  # Get path to reference pickle file.
  def refPicklPath(self):
    return IO.join(self.folder_name,"ref.pickl")

  # See if already scored.
  def isScored(self):
    return IO.isfile(self.refPicklPath())

  # Open reference pickle file and retrieve score information.
  def getScores(self,onlyAvg=False):
    reffile = open(self.refPicklPath())
    scores, avgscr, highpv, rank, sctype = load(reffile)
    reffile.close()
    if onlyAvg: return avgscr
    mincol,maxcol = -1,-1
    flatten = [scores[x] for x in scores if scores[x] != None]
    if sctype != 'RMSD' and sctype != 'RRMSD':
      flatten = [(1-x) for x in flatten] # Flip directions of floats.
    if len(flatten) != 0:
      minval,maxval = min(flatten),max(flatten)
      for key, val in scores.iteritems():
        if val == minval:   mincol = key
        elif val == maxval: maxcol = key
    return (avgscr, mincol, maxcol, highpv, rank)


