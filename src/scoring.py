''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++

scoring.py
June 20th, 2014; Alex Safatli

Handles Python wrapping of scoring in ABeRMuSA. To be implemented into main software suite. '''

import math, cPickle
from os.path import isfile
from utils import PDBnet
from utils import homology
from scipy.stats.kde import gaussian_kde

SCORE_TYPES = ['RRMSD','TMscore','GDT','RMSD']

''' Perform scoring. '''

def normpdf(x,mean,stdv):
    
    ''' Calculate the p-value of a normal distribution. '''
    
    variance = float(stdv)**2
    denom = (2*math.pi*variance)**(0.5)
    pv = (math.exp(-(float(x)-float(mean))**2/(2*variance)))/denom
    if pv < 1: return pv
    else: return 1.0


def score(aPDB,aFASTA,exe=None,logf=None):
    
    ''' Gets alignment score irregardless of alignment method. '''
    
    scores = []
    
    # Get PDB structure.
    p = PDBnet.PDBstructure(aPDB)
    
    # Get length of alignment.
    alignlen = len(p)
    
    # See what scores need to be done.
    if exe:
        scoresToDo = exe.scoresToDo
        if not scoresToDo: scoresToDo = SCORE_TYPES
    else: scoresToDo = SCORE_TYPES
    rrmsd,rpval,rmsd,tmsc,tpval,gdt = None,None,None,None,None,None
    
    # Get RRMSD and RMSD if length of alignment >= 100 residues.
    if 'RRMSD' in scoresToDo or 'RMSD' in scoresToDo:
        rrmsd, rmsd = homology.rrmsd(aPDB,aFASTA,True)
        if not exe.scpdbs or alignlen >= 100:
            rpval = normpdf(rrmsd,0.177,0.083)
        elif exe and logf and 'RRMSD' in scoresToDo:
            # Perform alignments in order to generate null distribution.
            logf.setTotalNum(logf.totalnum+2*(len(pdbli)+1))
            logf.writeTemporary(
                'Generating null distribution from SCOP for %s...' % (aPDB))
            scfolders = []
            run(exe.scpdbs,logf,ref=aPDB,exe=exe,quick=None)
            alignfldr = IO.getFileName(aPDB)
            o = open('%s/ref.pickl' % (alignfldr))
            dic, _, _ = cPickle.load(o)
            vals = dic.values()
            o.close()
            pdf = gaussian_kde(vals)
            rpval = pdf(rrmsd)
    
    # Get GDT and TMscore.
    if 'TMscore' in scoresToDo:
        tmsc = p.tmscore(aFASTA)
        tpval = 1 - math.exp(-math.exp((0.1512-tmsc)/0.0242))
         
    if 'GDT' in scoresToDo:
        gdt = p.gdt(aFASTA)
        
    # Add them to list in order as given.
    for it in scoresToDo:
        if it == 'RRMSD':
            scores.append(alignmentScore('RRMSD',rrmsd,rpval))
        elif it == 'RMSD':
            scores.append(alignmentScore('RMSD',rmsd))            
        elif it == 'TMscore':
            scores.append(alignmentScore('TMscore',tmsc,tpval))   
        elif it == 'GDT':
            scores.append(alignmentScore('GDT',gdt))
    
    # Return the scoring values.
    return scores

''' Score components for testing goodness of an alignment. '''

class alignmentScore(object):
    
    name,score,pval = '',0.0,None
    
    def __init__(self,key,val,pval=None):
        
        self.name  = key
        self.score = val
        self.pval  = pval
        
    def getScore(self):  return self.score
    def getName(self):   return self.name
    def getPValue(self): return self.pval
        
''' Score file. '''

class scoreFile(object):
    
    def __init__(self,path):
        
        self.path = path
        self.raw  = None
        self.scores = []
        if isfile(path):
            o = open(path,'r')
            self.raw = cPickle.load(o)
            o.close()
            self._read()
        
    def _read(self):
        
        if len(self.raw) == 3: # Legacy File
            sc,rc,pv = self.raw
            self.scores.append(alignmentScore('RRMSD',sc,pv))
            self.scores.append(alignmentScore('RMSD',rc))
        else:
            listofscores = self.raw
            for it in listofscores:
                key,val,pv = it
                self.scores.append(alignmentScore(key,val,pv))

    def iterScores(self):
        for it in self.scores:
            yield it.getName(), it.getScore(), it.getPValue()            
            
    def getScores(self):
        return [x for x in self.iterScores()]
    
    def getScoreByType(self,typ):
        for i in self.scores:
            if i.getName() == typ:
                return i.getName(), i.getScore(), i.getPValue()
    
    def addScores(self,items): self.scores.extend(items)
    
    def getRawReadFromFile(self): return self.raw
    
    def writeToFile(self):
        o = open(self.path,'w')
        scores = self.getScores()
        cPickle.dump(scores,o)
        o.close()
            