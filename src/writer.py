''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++

writer.py
Aug 28th, 2014; Alex Safatli

Construct a multiple sequence and structure alignment from all comprising pairwise alignments
that are done by ABeRMuSA. Creates a FASTA file with respective sequences for all PDBs and a
PDB file with models for each respective protein being aligned with Model 0 corresponding to the
reference protein. '''

import cPickle
from tempfile import NamedTemporaryFile as tempfile
from os.path import isfile
from os import unlink as delete
from shutil import move, copy
from numpy import mean
from subprocess import PIPE, Popen as system
from exewrapper import exeExists
from scoring import scoreFile
from utils import IO
from utils.homology import completePDB
from utils.logfile import XMLfile
from utils.FASTAnet import FASTAstructure as FASTA
from utils.PDBnet import PDBstructure as PDB

class profileAlignment(object):
    
    def __init__(self,args,startingFasta,keys,exe):
        
        ''' Given a pairwise FASTA sequence alignment, use MUSCLE to construct
        a profile alignment in a progressive manner. '''
        
        self.args     = args
        self.keys     = keys
        self.refname  = ''
        self.fastas   = [x for x in args if x != startingFasta]
        self.starting = startingFasta
        self.current  = None
        self.exe      = exe
        if not exeExists('muscle'):
            raise EnvironmentError('Executable muscle not present on system.')
        self._getReferenceName()
        
    def _getReferenceName(self):
        
        # Acquire the reference name from one of the arguments.
        
        fasta = FASTA(self.args[0],uniqueOnly=False)
        ref_pos = self.exe.plugin.ref_pos
        self.refname = fasta.orderedSequences[ref_pos].name        
    
    def _pairwise(self,f): 
        
        # Remove the reference from this FASTA if necessary.

        fasta = FASTA(f,uniqueOnly=False)
        fasta.removeSequence(self.refname)
        temp = tempfile(delete=False)
        fasta.writeFile(temp.name)
        f = temp.name
        if self.current != None: 
            compare = self.current
            outf    = self.current
        else:
            compare = self.starting
            temp = tempfile(delete=False)
            outf = temp.name
        system('muscle -profile -in1 %s -in2 %s -out %s' % (compare,f,outf),
               stdout=PIPE,stderr=PIPE,shell=True).communicate()[0]
        delete(f)
        if self.current is None: self.current = outf
        
    def write(self,fout):
        
        if len(self.fastas) == 0:
            # Only a single FASTA file to look at.
            copy(self.starting,fout)
        else:
            for fasta in self.fastas: self._pairwise(fasta)
            move(self.current,fout)
        self._postprocess(fout)
            
    def _postprocess(self,fout):
        
        # Open the resultant FASTA file and ensure all sequences are in position.
        fasta = FASTA(fout,uniqueOnly=False)
        fasta.reorderSequences([self.refname] + self.keys)
        fasta.writeFile(fout)

class multipleAlignment(object):

    def __init__(self,args,prefix,bestref,reffldr,logf,exe,
                 alphaC=False,optimize=False,touchup=False):
        
        self.args    = args # All comprising structures.
        self.prefix  = prefix # Label for output files.
        self.bestref = bestref # The reference to align all to.
        self.reffldr = reffldr # The folder for that reference.
        self.logf    = logf # The logfile object.
        self.scores  = {} # All scores associated with structures.
        self.exe     = exe # The exewrapper object.
        self.alphaC  = alphaC # Whether or not to use alpha-carbons.
        self.optim   = optimize # Whether or not to optimize.
    
    def _logParameters(self):
        
        ''' PRIVATE. Report on all input parameters. '''
        
        if self.alphaC:
            self.logf.write('NOTE; GM will contain Alpha Carbon coordinates.')
        if self.reffldr: self.logf.write('Reference folder <%s>' % (self.reffldr))
        if self.bestref: self.logf.write('Structure file <%s>' % (self.bestref[0]))
    
    def _getFASTAFiles(self):
        
        ''' PRIVATE. Acquire all FASTA files. '''
        
        fastas = []
        keys   = []
        
        # Determine the extension for the FASTA file (plugin-dependent).
        if self.exe:
            f_e = self.exe.plugin.fasta_ext
            cmd = self.exe.plugin.default_exe
        else: f_e = 'fasta'
        
        # See if arguments were provided explicitly.
        if self.args and self.exe:  
            # Argument list was provided.
            for f in self.args:
                key = IO.getFileName(f)
                # Is this the best reference?
                if key == self.reffldr: continue # Not aligned with itself.
                fasta = '%s/%s-%s.%s.%s' % (self.reffldr,cmd,self.reffldr,key,f_e)
                if isfile(fasta):
                    fastas.append(fasta)
                    keys.append(key)
                else: self.logf.write('WARNING; Missing file %s not included in GM.' % (fasta))
        else:
            # No argument list provided.
            manif = '%s/manifest.xml' % (self.reffldr)
            if not isfile(manif):
                # No XML manifest file. Incredibly old ABeRMuSA output.
                self.logf.write('WARNING; No manifest file was present. Possibly old ABeRMuSA output.')
                self.logf.write('WARNING; Legacy support of reference pickle files being resorted to.')
                if not isfile('%s/ref.pickl' % (self.reffldr)):
                    self.logf.write('ERROR; No reference pickle file present!')
                    return None, None
                o = open('%s/ref.pickl' % (self.reffldr))
                scores,_,_,_,_ = cPickle.load(o)
                o.close()
                for it in scores:
                    fasta = '%s/%s.%s' % (self.reffldr,it,f_e)
                    if isfile(fasta):
                        fastas.append(fasta)
                        keys.append(key)
                    else: self.logf.write('WARNING; Missing file %s not included in GM.' % (fasta))
            else:
                # An XML manifest file was present. Can acquire information from this.
                xml = XMLfile(manif)
                xml.read()
                if xml.root.tag != 'reference':
                    self.logf.write('ERROR; Invalid manifest file %s.' % (manif))
                    return None, None
                for it in xml.root:
                    if it.root == 'succeeded':
                        fasta = '%s/%s.%s' % (self.reffldr,it.text,f_e)
                        if isfile(fasta):
                            fastas.append(fasta)
                            keys.append(key)
                        else: self.logf.write('WARNING; Missing file %s not included in GM.' % (fasta))
                        
        return fastas, keys
                    
    def _filterNonSignificant(self,files,keys):
        
        bad = []
        for fi in files:
            name = IO.getFileName(fi)
            picklf = '%s/%s.pickl' % (self.reffldr,name)
            if not isfile(picklf):
                bad.append(fi)
                self.logf.write('Ignoring because not scored: <%s>.' % (name))
                continue
            scf = scoreFile(picklf)
            self.scores[fi] = scf.getScores()
            rr = scf.getScoreByType('RRMSD')
            if rr: sn,sc,pv = rr
            else: continue
            if pv < 0.05:
                bad.append(fi)
                self.logf.write('Omitting for P-value %f: <%s>.' % (pv,name))
         
        filtered = []
        keysout  = []
        for x in xrange(len(files)):
            if not files[x] in bad:
                filtered.append(files[x])
                keysout.append(keys[x])
        return filtered, keysout
    
    def _getLandmarks(self,files):
        
        ref_pos = self.exe.plugin.ref_pos
        ldata = {} # landmark lists
        ldata[self.reffldr] = []
        alnlen = -1        
        for fi in files:
            seqs = FASTAnet.FASTAstructure(fi,uniqueOnly=False)
            if len(seqs.sequences) != 2:
                self.logf.write('WARNING; <%s> does not contain exactly 2 sequences.' % (fi))
            if len(seqs.sequences) < 2: continue
            refstruc  = seqs.orderedSequences[ref_pos].sequence
            teststruc = seqs.orderedSequences[(ref_pos+1)%2].sequence
            if (teststruc == refstruc):
                self.logf.write('NOTE; <%s> contains identical sequences.' % (fi))
            elif len(refstruc) != len(teststruc):
                self.logf.write('ERROR; <%s> has non-matching lengths.' % (fi))
                return None
            f = IO.getFileName(fi) # Get name of PDB.
            ldata[f] = []
            if not ldata[self.reffldr]:
                for i in range(len(refstruc)):
                    if refstruc[i] != '-': ldata[self.reffldr].append(1) 
            for i in range(len(refstruc)):
                if refstruc[i] != '-':
                    if teststruc[i] == '-': ldata[f].append(0)
                    else: ldata[f].append(1)
            if alnlen == -1: alnlen = len(ldata[f])            
    
        return ldata    
    
    def _optimize(self,files,keys):
        
        self.logf.write('Optimizing for landmarks...')
        ldata = self._getLandmarks(files)
        if len(ldata.keys())-1 <= 0:
            self.logf.write('ERROR; Landmark list too small. Skipping optimization.')
            return files
        optli = {}
        for key in ldata.keys():
            if key == best: continue
            else: optli[key] = ldata[key]
        t = treeoptimizer(optli,bad,n-1,self.logf)        
        totalNum = t.getNumlandmarks()
        self.logf.write('Number of landmarks in total before: %d.' % (totalNum))
        t.optimize()
        totalNum = t.getNumLandmarks()
        o_bad = t.getOmitList()        
        for bad in o_bad:
            self.logf.write('Omitting for bad landmarks: <%s>.' % (bad))
        self.logf.write('Number of landmarks in total after: %d.' % (totalNum))
        
        filtered = []
        keysout  = []
        for x in xrange(len(files)):
            if not IO.getFileName(files[x]) in o_bad:
                filtered.append(files[x])
                keysout.append(keys[x])
        return filtered, keysout
    
    def _checkPDB(self,pdb):
        
        if pdb.CheckComplete() != True:
            self.logf.write('Found incomplete PDB data in %s. Completing.' % (
                pdb.filepath))
            temp = tempfile(delete=False)
            pdb.write(temp.name)
            completePDB(temp.name,pdb.filepath)
            delete(temp.name)
            pdb = PDB(pdb.filepath)
        return pdb
    
    def _getOriginalIndices(self,key):
        
        pdb = ''
        for arg in self.args:
            if key in arg:
                pdb = arg
                break
        if pdb != '' and isfile(pdb):
            p = PDB(pdb)
            frstch = p.chains.keys()[0]
            return p.GetChain(frstch).GetIndices()
        else: return None
    
    def _integrateChainToAlignedPDB(self,p,o,ch,index,key):
        
        p.AddModelToStructure(index,o.GetChain(ch))
        residues = p.GetModel(index).GetResidues()
        originalindices = self._getOriginalIndices(key)
        numstart = 1
        for it in xrange(len(residues)):
            res = residues[it]
            res.chain = ''
            res.model = index
            res.index = originalindices[it]
            for at in res.GetAtoms():
                at.serial = str(numstart)
                numstart += 1
    
    def _preambleForPDB(self,files):
        
        liout = []
        st = 'Alignment automatically generated using ABeRMuSA:' + \
            ' Approximate Best  Reference Multiple Structure Alignment' + \
            ' ver. %s.' % (self.exe.version)
        liout.append(st)
        liout.append('')
        st = 'Reference Name: <%s> (MODEL 0)' % (self.reffldr)
        liout.append(st)
        scores = {}
        for fi in files:
            for sn,sc,pv in self.scores[fi]:
                if not sn in scores: scores[sn] = []
                scores[sn].append(sc)
        for sn in scores:
            liout.append('Average Score <%s> = %.3f' % (sn,mean(scores[sn])))
        liout.append('')
        liout.append('File/Chain/Model No. Mappings')
        return liout    
    
    def construct(self):
        
        ''' Construct the alignment. '''
        
        # Report and get FASTA files.
        self._logParameters()
        fastas, keys = self._getFASTAFiles()
        if fastas == None: return None
        fastas, keys = self._filterNonSignificant(fastas, keys)
        if self.optim: fastas, keys = self._optimize(fastas)
        
        # Make a single sequence alignment.
        self.logf.write('Consolidating all FASTA files into a single alignment...')
        multi = profileAlignment(fastas,fastas[0],keys,self.exe)
        multi.write('%s.fasta' % (self.prefix))
        self.logf.write('Wrote multiple sequence alignment of structures to %s.' % (
            self.prefix + '.fasta'))
        
        # Make a single PDB file.
        p = PDB()
        for re in self._preambleForPDB(fastas): p.AddRemark(re)
        m = open('%s.pdb_map' % (self.prefix),'w')
        self.logf.write('Consolidating all PDB data into a single PDB file...')
        chs,ref_pos,index = ['A','B'],self.exe.plugin.ref_pos,0
        for x in xrange(len(fastas)):
            fi = fastas[x]
            key = keys[x]
            nm = '%s/%s.pdb' % (self.reffldr,IO.getFileName(fi))
            pdbdata = PDB(nm)
            pdbdata = self._checkPDB(pdbdata)
            if index == 0:
                # Add reference only once.
                self._integrateChainToAlignedPDB(p,pdbdata,chs[ref_pos],index,key)
                m.write('%s\t%s\tMODEL %d\n' % (nm,chs[ref_pos],index))
                p.AddRemark('%s  %s  MODEL %d' % (nm,chs[ref_pos],index))
                index += 1
            self._integrateChainToAlignedPDB(p,pdbdata,chs[(ref_pos+1)%2],index,key)
            m.write('%s\t%s\tMODEL %d\n' % (nm,chs[(ref_pos+1)%2],index))
            p.AddRemark('%s  %s  MODEL %d' % (nm,chs[(ref_pos+1)%2],index))
            index += 1
        p.write('%s.pdb' % (self.prefix))
        m.close()
        self.logf.write('Wrote multiple structure alignment of structures to %s.' % (
            self.prefix + '.pdb'))
        
        # Write GM file and landmark file.
        p.WriteGM('%s.fasta' % (self.prefix),'%s.gm' % (
            self.prefix),ismodel=True,CA=self.alphaC)
        self.logf.write('Wrote GM data successfully to %s.' % ('%s.gm' % (self.prefix)))
        p.WriteLandmarks('%s.fasta' % (self.prefix),'%s.landmarks' % (self.prefix),
                         ismodel=True)
        self.logf.write('Wrote landmark file successfully to %s.' % (
            '%s.landmarks' % (self.prefix)))
        
        return ['%s.gm' % (self.prefix),'%s.landmarks' % (self.prefix),
                '%s.fasta' % (self.prefix),'%s.pdb' % (self.prefix)]