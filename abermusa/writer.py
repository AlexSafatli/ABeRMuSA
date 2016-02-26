''' Construct a multiple sequence and structure alignment from all comprising pairwise alignments
that are done by ABeRMuSA. Creates a FASTA file with respective sequences for all PDBs and a
PDB file with models for each respective protein being aligned with Model 0 corresponding to the
reference protein. '''

# Date:   Aug 28 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import cPickle
import numpy as np
from tempfile import NamedTemporaryFile as tempfile
from os.path import isfile
from os import unlink as delete
from shutil import move, copy
from numpy import mean
from subprocess import PIPE, Popen as system
from exewrapper import exeExists
from scoring import scoreFile
from labblouin import IO
from labblouin.homology import completePDB
from labblouin.logfile import XMLfile, logfile
from labblouin.FASTAnet import FASTAstructure as FASTA
from labblouin.PDBnet import PDBstructure as PDB
from treeoptimizer import treeoptimizer


class profileAlignment(object):

	def __init__(self,args,startingFasta,keys,reffldr,exe):

		''' Given a pairwise FASTA sequence alignment, use MUSCLE to construct
        a profile alignment in a progressive manner. '''

		self.args     = args
		self.keys     = keys
		self.fastas   = [x for x in args if x != startingFasta]
		self.refname  = reffldr
		self.exe      = exe
		self.handles  = self._getHandles()
		self.starting = startingFasta
		self.current  = None
		if not exeExists('muscle'):
			raise EnvironmentError('Executable "muscle" was not found on path.')

	def getFileKey(self,fipath):

		# Given a filepath, return the key corresponding to it.

		return self.keys[self.args.index(fipath)]

	def getHandle(self,fipath):

		# Given a filepath, return the FASTA handle.

		return self.handles[self.getFileKey(fipath)]

	def _getHandles(self):

		# Get all FASTAnet handles.

		dicout  = {}
		ref_pos = self.exe.plugin.ref_pos
		for f in xrange(len(self.args)):
			# Get the handle.
			dicout[self.keys[f]] = FASTA(self.args[f],uniqueOnly=False)
			# Correct the names in that handle to match keys.
			oldname = lambda d: dicout[self.keys[f]].orderedSequences[d].name
			dicout[self.keys[f]].renameSequence(oldname(ref_pos),self.refname)
			dicout[self.keys[f]].renameSequence(oldname((ref_pos+1)%2),self.keys[f])
		return dicout

	def _pairwise(self,f,refine=True): 

		# Remove the reference from this FASTA if necessary.
		fasta = self.getHandle(f)
		fasta.removeSequence(self.refname)
		temp = tempfile(delete=False)
		fasta.writeFile(temp.name)
		f = temp.name
		if self.current != None: 
			compare = self.current
			outf    = self.current
		else:
			start = self.getHandle(self.starting)
			stemp = tempfile(delete=False)
			start.writeFile(stemp.name)
			compare = stemp.name
			temp = tempfile(delete=False)
			outf = temp.name
		system('muscle -profile -in1 %s -in2 %s -out %s' % (compare,f,outf),
		       stdout=PIPE,stderr=PIPE,shell=True).communicate()[0]
		delete(f)
		# refine the profile
		if refine:
			system('muscle -in %s -out %s -refine'% (outf,outf),
			       stdout=PIPE,stderr=PIPE,shell=True).communicate()[0]
		if self.current is None:
			delete(stemp.name)
			self.current = outf


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


class Pairwise2Multiple:
	''' Writes a GM file based on pairwise alignments'''

	def __init__(self,prefix,fastas,alnpdb,algnr,log,atomtype='centroid'):
		self.log = log
		self.aligner = algnr
		self.fastas=fastas
		self.atomtype=atomtype
		if isinstance(alnpdb,str): self.pdb = PDB(alnpdb)
		elif isinstance(alnpdb,PDB): self.pdb = alnpdb
		self.prefix = prefix
		self.parseRemarks()
		self.PairwiseHomologs()
		self.GMconglomerate()
		np.savetxt(self.prefix+'.gm', self.gm, fmt='%s', delimiter=';')
		self.log.write('GM file written to %s, with %d rows and %d columns'%(self.prefix+".gm",self.gm.shape[0],self.gm.shape[1]-1))

	def parseRemarks(self):
		'''get equivalences from an ABeRMuSA aligned pdb'''
		eq={}
		eqfiles={}
		mapfile = self.prefix+'.pdb_map'
		ref=None
		if isfile(mapfile):
			with open(mapfile) as mp:
				for line in mp:
					bl=line.strip().split()
					fil = bl[0].strip('.pdb')
					if int(bl[-1]) == 0:
						ref = bl[0].split('/')[0]
						eq[0] = ref
						eqfiles[fil]=0
					else:
						name  = fil.split('/%s-'%(self.aligner))[-1]
						if ref:
							name = name.split('.')[-1]
						model = int(bl[-1])
						eq[model]=name
						eqfiles[fil]=model
		else:
			rmks=self.pdb.GetRemarks()
			ref=None
			for i in rmks:
				if "Reference Name:" in i:
					ref = i[i.find('<')+1:i.find('>')].strip('.pdb').strip('.fasta')
					eq[0]=ref
				elif (ref in i) and ("Reference Name:" not in i):
					bl = i.strip().split()
					if ' '.join(bl[-2:]) == 'MODEL 0': continue
					model = int(bl[-1])
					eqfiles[bl[1].strip('.pdb').strip('fasta')]=model
					if ref:
						name = bl[1].strip('%s/%s-%s'%(ref,self.aligner,ref)).strip('.pdb').strip('.fasta')
					#name = bl[1].split('.')[-2]
					eq[model]=name
		self.eqpdb=eq
		self.eqfiles=eqfiles


	def PairwiseHomologs(self):
		''' 
		Giving a list of fasta files return a dictionary with a list of tuples of aligned indices.
		the dicionary key will be the fasta filename
		'''
		fastas = self.fastas
		d={}
		for fasta in fastas:
			fasta = fasta.strip('.fasta').strip('.pdb')
			bl = fasta.split('/%s-'%(self.aligner))
			refname = bl[0]
			tarname = bl[1]
			fas = FASTA(fasta+'.fasta',uniqueOnly=False)
			hom = fas.findAlignedResidueIndices()
			seqs= fas.sequences
			#ref=[]
			#tar=[]
			for s in seqs:
				if refname in s:
					ref = seqs[s].toIndices()
					ref = np.array(ref)[hom].astype('int').tolist()
					#ref.append(re)
				else:
					tar = seqs[s].toIndices()
					tar = np.array(tar)[hom].astype('int').tolist()
					#tar.append(ta)
			#d[fasta]=zip(ref[-1],tar[-1])
			d[fasta]=zip(ref,tar)
		self.hom = d
		self.refname = refname


	def GMconglomerate(self):
		'''
		giving the intersection list of fully homologous reference, and the equivalence tuples
		Generate a gm file and a landmarks file
		'''
		self.conglomerateRefhomologs()
		refmod = self.pdb.GetModel(0)
		refgm=[]
		names=[]
		land={self.refname:[]}
		for r in self.cong:
			res = refmod.GetResidueByPosition(r)
			if self.atomtype =='centroid': cent= res.Centroid()
			else: cent = res.GetCA()
			refgm.extend([cent.x,cent.y,cent.z])
			land[self.refname].append((str(self.cong.index(r)),str(res.index),str(res.name)))
		gm = np.array(refgm)
		names.append('>%s'%(self.refname))
		count = 0
		for h in sorted(self.hom):
			h = h.strip('.pdb').strip('.fasta')
			if h in self.eqfiles:
				temp=[]
				m = self.pdb.GetModel(self.eqfiles[h])
				n = h.split('.')[-1]
				land[n]=[]
				index = -1
				for ind in self.hom[h]:
					if ind[0] in self.cong:
						index +=1
						res = m.GetResidueByPosition(ind[1])
						if self.atomtype =='centroid': cent= res.Centroid()
						else: cent = res.GetCA()
						temp.extend([cent.x,cent.y,cent.z])
						land[n].append((str(index),str(res.index),str(res.name)))
				names.append('>'+self.eqpdb[self.eqfiles[h]])
				gm = np.vstack((gm,np.array(temp)))
				count+=1
		#print 'number of processed fasta: %d'%(count)
		# write landmark file
		line=''
		for k,v in land.iteritems():
			line+='>%s\n'%(k)
			line+= '\n'.join(['\t'.join(x) for x in v])+'\n'
		with open(self.prefix+'.landmarks','w') as F: F.write(line[:-1])
		self.gm = np.column_stack((np.array(names), gm))
		#print 'Dimensions of gm matrix: %d, %d'%(self.gm.shape)

	def conglomerateRefhomologs(self):
		''' create an intersection list of fully homologous reference'''
		hom = self.hom
		cong=[x[0] for x in hom[hom.keys()[0]]]
		for k,v in hom.iteritems():
			a = [x[0] for x in v]
			cong=list(set(cong).intersection(a))
		self.cong=cong

class multipleAlignment(object):

	def __init__(self,args,prefix,bestref,reffldr,logf,exe,
	             atomtype='centroid',curate=False,optimize=False,touchup=False,
	             MD=False,fasta=True,gm=False):

		self.args    = args # All comprising structures.
		self.prefix  = prefix # Label for output files.
		self.bestref = bestref # The reference to align all to.
		self.reffldr = reffldr # The folder for that reference.
		self.logf    = logf # The logfile object.
		self.scores  = {} # All scores associated with structures.
		self.exe     = exe # The exewrapper object.
		self.atomtype  = atomtype # Whether or not to use alpha-carbons.
		self.curate  = curate # Whether or not to curate output PDB files.
		self.optim   = optimize # Whether or not to optimize.
		self.MD      = MD # Whether or not to assume a trajectory from MD.
		self.fasta   = fasta # Whether or not to write a FASTA file.
		self.gm      = gm # Whether or not to write a gm file NOT based on fasta.

	def _logParameters(self):

		''' PRIVATE. Report on all input parameters. '''

		if self.atomtype == 'CA':
			self.logf.write('Alpha Carbon coordinates in GM (atomtype="CA").', status='NOTE')
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
				else: self.logf.write('Missing file %s not included in GM.' %
				                      (fasta), status='WARNING')
		else:
			# No argument list provided.
			manif = '%s/manifest.xml' % (self.reffldr)
			if not isfile(manif):
				# No XML manifest file. Incredibly old ABeRMuSA output.
				self.logf.write('No manifest file (%s) was present!'%(manif), status='WARNING')
				self.logf.write('Using legacy support of ref. pickle files.', status='NOTE')
				refpicklfi = '%s/ref.pickl' % (self.reffldr)
				if not isfile(refpicklfi):
					self.logf.write('No reference pickle file (%s) present!' % (
					    refpicklfi), status='ERROR')
					return None, None
				o = open('%s/ref.pickl' % (self.reffldr))
				scores,_,_,_,_ = cPickle.load(o)
				o.close()
				for it in scores:
					fasta = '%s/%s.%s' % (self.reffldr,it,f_e)
					if isfile(fasta):
						fastas.append(fasta)
						keys.append(key)
					else: self.logf.write('Missing file %s not incl. in GM.' 
					                      % (fasta), status='ERROR')
			else:
				# An XML manifest file was present. Can acquire information from this.
				xml = XMLfile(manif)
				xml.read()
				if xml.root.tag != 'reference':
					self.logf.write('Invalid manifest file %s.' % (manif), status='ERROR')
					return None, None
				for it in xml.root:
					if it.root == 'succeeded':
						fasta = '%s/%s.%s' % (self.reffldr,it.text,f_e)
						if isfile(fasta):
							fastas.append(fasta)
							keys.append(key)
						else: self.logf.write('Missing file %s not incl. in GM.' 
						                      % (fasta), status='WARNING')

		return fastas, keys

	def _filterNonSignificant(self,files,keys):

		self.logf.write('Filtering non-significant files (len=%d) from alignment ' 
		                % (len(files)) + '(using RRMSD)...')
		bad = []
		for fi in files:
			name = IO.getFileName(fi)
			picklf = '%s/%s.pickl' % (self.reffldr,name)
			if not isfile(picklf):
				bad.append(fi)
				self.logf.write('Ignoring because no pickle file (%s) found: <%s>.' %(
				    picklf,name), status='NOTE')
				continue
			scf = scoreFile(picklf)
			self.scores[fi] = scf.getScores()
			rr = scf.getScoreByType('RRMSD')
			if rr: sn,sc,pv = rr
			else: continue
			if pv > 0.95:
				bad.append(fi)
				self.logf.write('Omitting for P-value %f: <%s>.' % (pv,name))

		filtered = []
		keysout  = []
		for x in xrange(len(files)):
			if not files[x] in bad:
				filtered.append(files[x])
				keysout.append(keys[x])
		self.logf.write('Filtered list of files has %d FASTA files.' % (
		    len(filtered)))
		return filtered, keysout, bad

	def _getLandmarks(self,files):

		ref_pos = self.exe.plugin.ref_pos
		ldata = {} # landmark lists
		ldata[self.reffldr] = []
		alnlen = -1        
		for fi in files:
			seqs = FASTA(fi,uniqueOnly=False)
			if len(seqs.sequences) != 2:
				self.logf.write('<%s> does not contain exactly 2 sequences.' % (fi),
					status='WARNING')
			if len(seqs.sequences) < 2: continue
			refstruc  = seqs.orderedSequences[ref_pos].sequence
			teststruc = seqs.orderedSequences[(ref_pos+1)%2].sequence
			if (teststruc == refstruc):
				self.logf.write('<%s> contains identical sequences.' % (fi), status='NOTE')
			elif len(refstruc) != len(teststruc):
				self.logf.write('<%s> has non-matching lengths.' % (fi), status='ERROR')
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

	def _optimize(self,files,keys,bad):

		self.logf.write('Optimizing for landmarks...')
		ldata = self._getLandmarks(files)
		if len(ldata.keys())-1 <= 0:
			self.logf.write('Landmark list too small. Skipping optimization.', status='ERROR')
			return files
		optli = {}
		for key in ldata.keys():
			if key == self.bestref: continue#if key == best: continue
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
			    pdb.filepath), status='NOTE')
			temp = tempfile(delete=False)
			pdb.WriteFile(temp.name)
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
		else:
			raise IOError('Could not find original PDB for %s amongst arguments.' % 
			              (key))

	def _integrateChainToAlignedPDB(self,final,temp,o,ch,index,key):

		# Add model to temporary structure.
		temp.AddModel(index,o.GetChain(ch).AsModel(index))

		# Fix all residues and indices.
		residues = temp.GetModel(index).GetResidues()
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

		# Write this model to the file (in append mode).
		o = open(final,'a')
		o.write(str(temp.GetModel(index)))
		o.close()

		# Remove changes made to temporary structure.
		temp.RemoveModel(index)

	def _preambleForPDB(self,files):

		liout = []
		st = 'Alignment automatically generated using ABeRMuSA:' + \
		    ' Approximate Best Reference Multiple Structure Alignment' + \
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

	def _alignment4MD(self,fastas):

		''' Constructs a dummy alignment file when input is a trajectory. '''

		fout= open('%s.fasta'%(self.prefix),'w')
		seq = '\n'.join(open(fastas[0]).read().split('>')[1].split()[1:])
		for i in range(len(fastas)+1): # All the models plus the "reference".
			fout.write('>Model%d\n%s\n'%(i,seq))

	def construct(self):

		''' Construct the alignment. '''

		# Track what files are written.
		filesWritten = []

		# Report and get FASTA files.
		self._logParameters()
		fastas, keys = self._getFASTAFiles()
		if fastas == None: return None
		fastas, keys, bad = self._filterNonSignificant(fastas, keys)
		if self.optim: fastas, keys = self._optimize(fastas,keys,bad)

		# Make a single sequence alignment.
		if not self.MD and self.fasta:
			self.logf.write('Consolidating all FASTA files into an alignment...')
			if isfile(self.prefix+'.fasta'): self.logf.write('Fasta file already in path')
			else:
				multi = profileAlignment(fastas,fastas[0],keys,self.reffldr,self.exe)
				multi.write('%s.fasta' % (self.prefix))
				filesWritten.append('%s.fasta' % (self.prefix))
				self.logf.write('Wrote multiple sequence alignment FASTA to %s.' % (
				    self.prefix + '.fasta'))

		elif self.MD and self.fasta:
			self._alignment4MD(fastas)
			self.logf.write('Wrote dummy multiple sequence alignment FASTA (MD) to %s.' % (
			    self.prefix + '.fasta')) 
			filesWritten.append('%s.fasta' % (self.prefix))

		elif not self.fasta and not self.gm:
			self.logf.write('Skipping the writing of MSA to FASTA.')

		# Report on writing status.
		self.logf.write('Consolidating all PDB data into a single PDB file...')

		# Acquire preamble for remarks.
		finalRemarks = [re for re in self._preambleForPDB(fastas)]

		# Construct a PDB map file.
		m = open('%s.pdb_map' % (self.prefix),'w')
		filesWritten.append('%s.pdb_map' % (self.prefix))

		# Make a single PDB file.
		temporaryPDB = PDB() # Empty PDB structure for manipulation.
		tempPath = '%s_part.pdb' % (self.prefix)
		finalPath = '%s.pdb' % (self.prefix)
		chs,ref_pos,index = ['A','B'],self.exe.plugin.ref_pos,0
		filesWritten.append(finalPath)

		# Read all FASTA files.
		numFastas = len(fastas)
		self.logf.write('Reading %d FASTA files (after filtering)...'%(numFastas))
		for x,fi in enumerate(fastas):
		#for fi,x in enumerate(fastas):
			# Get file data.
			key = keys[x]
			nm  = '%s/%s.pdb' % (self.reffldr,IO.getFileName(fi))
			pdbdata = PDB(nm)
			if self.curate: pdbdata = self._checkPDB(pdbdata)

			if index == 0:
				# Add reference only once.
				self._integrateChainToAlignedPDB(tempPath,temporaryPDB,pdbdata,
				                                 chs[ref_pos],index,self.reffldr)
				m.write('%s\t%s\tMODEL %d\n' % (nm,chs[ref_pos],index))
				finalRemarks.append('%s  %s  MODEL %d' % (nm,chs[ref_pos],index))
				index += 1

			self._integrateChainToAlignedPDB(tempPath,temporaryPDB,pdbdata,
			                                 chs[(ref_pos+1)%2],index,key)
			m.write('%s\t%s\tMODEL %d\n' % (nm,chs[(ref_pos+1)%2],index))
			finalRemarks.append('%s  %s  MODEL %d' % (nm,chs[(ref_pos+1)%2],index))
			index += 1

		m.close()

		remarkStr = ''
		for it in xrange(len(finalRemarks)):
			remarkStr += 'REMARK%4s %s\n' % (str(it+1),finalRemarks[it])        
		q = open(finalPath,'w')
		r = open(tempPath,'r')
		q.write(remarkStr)
		q.write(r.read())
		q.close()
		r.close()            
		delete(tempPath)
		self.logf.write('Wrote multiple structure alignment of structures to %s.' % (
		    self.prefix + '.pdb'))


		# Write GM file and landmark file.
		if self.fasta and not gm:
			p = PDB(finalPath)
			p.WriteGM('%s.fasta' % (self.prefix),'%s.gm' % (
			    self.prefix),atomtype=self.atomtype)
			self.logf.write('Wrote GM data successfully to %s.' % ('%s.gm' % (self.prefix)))
			p.WriteLandmarks('%s.fasta' % (self.prefix),'%s.landmarks' % (self.prefix))
			self.logf.write('Wrote landmark file successfully to %s.' % (
			    '%s.landmarks' % (self.prefix)))
			filesWritten.append('%s.gm' % (self.prefix))
			filesWritten.append('%s.landmarks' % (self.prefix))

		elif not self.fasta and self.gm:
			algnr = self.exe.plugin.default_exe
			alnpdb = PDB(finalPath)
			multi = Pairwise2Multiple(self.prefix, fastas, alnpdb, algnr, self.logf,atomtype=self.atomtype)
			self.logf.write('Skipping the writing of MSA to FASTA, but writing GM based on pairwise alignments.')
			filesWritten.append('%s.gm' % (self.prefix))
			filesWritten.append('%s.landmarks' % (self.prefix))			

		else: self.logf.write('Skipping writing of GM and landmark file because FASTA was not written, and GM was not requested.')

		return filesWritten

# For purposes of generating a multiple str. alignment RMSD.

def GMtoMatrix(prefix):
	f = open(prefix + '.gm')
	s = {} ## s for structures
	for line in f:
		if line.strip().endswith(';'):
			line = line.strip()[:-1]
		temp = [[], [], []]
		bline = line[line.find(';')+1:].strip().split(';')
		for i in range(len(bline)):
			temp[i%3].append(float(bline[i]))
		s[line[:line.find(';')]] = temp
	return s

def RMSDfromMEAN(p,m):
	# Calculates the RMSD for the meanshape (m) and the problem shape (p)
	sum_dist=0.0
	count = 0
	for j in range(len(m[0])):
		d =(m[0][j]-p[0][j])**2 + (m[1][j]-p[1][j])**2 + (m[2][j]-p[2][j])**2
		sum_dist += d
		count += 1
	# calculate the sqrt of the mean deviation
	RMSD = math.sqrt(sum_dist/count)

	return RMSD

class pseudoAlignment(multipleAlignment):

	def __init__(self,refldr):
		super(pseudoAlignment,self).__init__([],'temp',None,refldr,
		                                     logfile('',False,True),None)

	def getRMSD(self):
		matrix   = []
		outfiles = super(pseudoAlignment,self).construct()
		if outfiles: structures = GMtoMatrix(self.prefix)
		else:        return -1.0
		names = sorted(structures.keys())
		for i in names:
			for j in names:
				if i != j:
					matrix.append(RMSDfromMEAN(structures[j],structures[i]))
		for f in outfiles: os.remove(f) # Clean outfiles.
		return mean(matrix)
