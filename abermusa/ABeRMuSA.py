#!/usr/bin/env python

''' ABeRMuSA (Approximate Best Reference Multiple Structure Alignment)

Complete refactoring of code written by Khanh Nguyen in 2012.

Automatic pairwise alignment capable of being extended to any 
pairwise alignment executable by a plugin interface (see plugins/).

By default, this application uses MATT, a flexible pairwise aligner.

Input:   PDB files -OR- folders of PDB files (>= 2 PDBs).
Options: See help menu (--help, -h). 
Output:  A multiple structure alignment (PDB and FASTA file). '''

# Date: May 6 2013
# Main Author: Alex Safatli
# E-mail: safatli@cs.dal.ca
# Contributions: Jose Sergio Hleap, Khanh Nguyen

from __version__ import VERSION

__author__    = 'Alexander Safatli, Jose Sergio Hleap, Khanh Nguyen'
__email__     = 'safatli@cs.dal.ca, jshleap@squalus.org'
__version__   = VERSION

import numpy as np
import pandas as pd
from datetime import datetime
from shutil import rmtree as rmDir
from os import path, mkdir, rename, remove
from os.path import join, isfile, isdir
from shutil import copyfile, copy
from copy import deepcopy
import optparse, glob, sys, tarfile

import aligner
from labblouin import homology, pfam, scop, PDBnet, IO
from labblouin.logfile import logfile, XMLfile, timer as T
from exewrapper import exewrapper, exeExists as isCmd
from quickref import quickref
from writer import multipleAlignment
from scoring import SCORE_TYPES
from plugins import *

# Constants, Initialization

SCRIPT_FOLDER = path.split(path.realpath(__file__))[0]
PLUGIN_FOLDER = path.join(SCRIPT_FOLDER,'plugins')
PLUGIN_PYS    = glob.glob(path.join(PLUGIN_FOLDER,'*.py'))
PLUGINS       = [path.split(x)[-1].strip('.py') for x in PLUGIN_PYS \
                 if not x.endswith('__init__.py')]
PDB_ALLOW     = ['pdb','ent','atm'] # Extensions to consider for PDB files.
PDB_FOLDER    = '_input' # Folder to create to store input PDB files.
PDB_CACHE     = '_scop' # Folder to create to store files downloaded from SCOP.

# Auxiliary Functions, Classes

class tarArchive:
    def __init__(self,prefix):
        self.prefix = prefix
        self.fname  = prefix + '.tar.gz'
        self.tar    = tarfile.open(self.fname,"w:gz")
    def add(self,fi):
        self.tar.add(fi)
        return fi
    def close(self):
        self.tar.close()

class referenceCollection(object):
    def __init__(self,args):
        self.references = []
        self.foundli    = []
        for it in args:
            self.references.append(it)
            self.foundli.append(False)
    def get(self): return self.references
    def isFound(self,ref):
        return self.foundli[self.references.index(ref)]
    def found(self,ref):
        self.foundli[self.references.index(ref)] = True
    def map(self,func): map(func, self.references)
    def setReference(self,ref,newref):
        self.references[self.references.index(ref)] = newref
    def hasReference(self): return (len(self.references) > 0)
    def iterNotFound(self):
        for i in xrange(len(self.foundli)):
            if self.foundli[i] == None:
                yield self.references[i]    
    def anyNotFound(self):
        for found in self.foundli:
            if found == False:
                return True
        return False

def stripAtomsFromPDB(fi,fn):

    ''' Get atomlines from fi file and write to fn. (JSH). '''

    if path.isfile(fn): return fn
    o = open(fn,'w')
    with open(fi) as F:
        for line in F:
            if line.startswith('ATOM'): o.write(line)
        o.close()
    return fn

def enumerateFilesInFolder(path):

    ''' Given a directory, look in that directory for any
    alignable items. '''

    files = []

    # If there is a reference pickle file in path already,
    # already aligned.
    refpi = join(path,'ref.pickl')
    if not isfile(refpi):
        # Go over all allowable extensions.
        for ext in PDB_ALLOW:
            lookf = join(path,'*.' + ext)
            enum  = glob.glob(lookf)
            files.extend(enum)

    return files

def splitTraj(filename,log):

    ''' Split a trajectory PDB file (from a molecular dynamics
    simulation) into individual files (JSH). '''

    multi = []
    with open(filename) as F:
        model = ''
        modelno = 0
        modelfi = ''
        for line in F:
            if line.startswith('MODEL'):
                bl = line.strip().split()
                modelno = int(bl[1])
                modelfi = path.join(PDB_FOLDER,'Model%d.pdb'%(modelno))
            elif line.startswith('ATOM'): model+=line
            elif line.startswith('ENDMDL'):
                if model != '':                           
                    if not path.isfile(modelfi):
                        log.writeTemporary('Extracted <%s> (Model %s) to <%s>.' % (
                            filename,modelno,modelfi))
                        fout=open(modelfi,'w')
                        fout.write(model)
                        fout.close()
                        model = ''	
                        multi.append(modelfi)
                    else:
                        log.writeTemporary('File Model%s.pdb already in input folder.' % (modelno))
                        model = ''	
                        multi.append(modelfi)                
            else: continue 

    return multi

def handleFile(fi,log,refw,clean=False,split=False,MD=False):

    ''' Handle a file during startup; do cleaning, etc. '''

    # Is it the reference?
    isref, refkey = False, None
    if fi in refw.references:
        refw.found(fi)
        isref = True
        refkey = fi

    # Need only ATOM lines.
    fn = IO.getFileName(fi) # Get filename.
    if not MD:
        af = path.join(PDB_FOLDER,fn+'.pdb')
        fi = stripAtomsFromPDB(fi,af)
    if isref:
        refw.setReference(refkey,fi)
        refkey = fi

    # Translate to PDBnet structure if necessary.
    if (clean or split):
        try: p = PDBnet.PDBstructure(fi)
        except:
            log.write('<%s> could not be parsed as PDB format.' % (fn), status='WARNING')
            if isref:
                log.write('<%s> is provided reference and could not be parsed.' % (fn), status='ERROR')
                exit(1)
            return []

    # See if needs to be cleaned.
    if clean and not p.CheckComplete():
        cf = path.join(PDB_FOLDER,fn+'-c.pdb')
        if not path.isfile(cf): homology.completePDB(fi,cf)
        fi = cf
        if isref:
            refw.setReference(refkey,fi)    
            refkey = fi

    # Multiple chains?
    if split:
        chains = p.orderofchains
        numch  = len(chains)
        if (numch > 1):
            # Split it.
            log.writeTemporary('<%s> has multiple chains (%d). Extracting.' % (
                fn,numch), status='NOTE')
            multi = []
            for ch in chains:
                cf = path.join(PDB_FOLDER,fn+'_%s.pdb' %(ch))
                if not path.isfile(cf): pfam.extractPDBChain(fi,ch,cf)
                log.writeTemporary('Extracted <%s> (Chain %s) to <%s>.' % (fn,ch,cf))
                multi.append(cf)
            return multi

    # Multiple models?
    if MD:
        # Split it.
        log.write('NOTE <%s> has multiple models. Extracting.' % (fn))
        multi = splitTraj(fi,log)
        if refw.hasReference():
            refw.map(lambda t: refw.setReference(t,path.join(PDB_FOLDER,t)))
            refw.map(lambda t: refw.found(t))
        return multi

    # Return the file (as a single-element list).
    return [fi]

def acquireFiles(arg,fl,log,refw,clean=False,split=False,MD=False):

    ''' Populate an input filelist with alignable items from 
    the arguments. '''

    # See if any arguments at all.
    if len(arg) == 0: return # No arguments provided.

    # See if folder for PDBs has been made; otherwise, make it.
    if not path.isdir(PDB_FOLDER): mkdir(PDB_FOLDER)
    log.write('Created folder %s for post-processed (input) PDB files.' % (
        PDB_FOLDER))    

    # Enumerate over arguments.
    for fi in arg:
        # Is a directory?
        if path.isdir(fi):
            # Enumerate allowable items in directory.
            darg = enumerateFilesInFolder(fi)
            acquireFiles(darg,fl,log,refw,clean,split,MD)
        # Is a file?
        elif path.isfile(fi):
            fl.extend(handleFile(fi,log,refw,clean,split,MD))
        # Is not file BUT is reference?
        elif fi in refw.references:
            log.writeTemporary('<%s> set as reference but could not be found.' % (
                fi), status='ERROR')
            exit(1)
        # Not file or folder.
        else: log.writeTemporary('<%s> not a file or folder. Skipping argument.' % (fi), status='WARNING')

def Procrustes(prefix,traj,fasta=None,dim=3,atomtype='centroid',log=None):
    '''
    Perform partial procrustes superimposition on a pdb trajectory file (for MD only). Is based on the
    R package shapes (Dryden, 2014) through rpy2

    :param prefix: The prefix of outputs
    :type prefix: string
    :param traj: A Filename of trajectory of a MD simulation in PDB format
    :type traj: string
    '''
    line = 'If procrustes is used, please cite:\n'
    line+= '"Ian L. Dryden (2014). shapes: Statistical shape analysis. R package version 1.1-10."'
    line+= ' http://CRAN.R-project.org/package=shapes'
    if log: log.write(line,silent=False)
    else: print line
    prefix =prefix+'_gps'
    from labblouin.rpy2GPA import Rmatrix2Rarray,array2R,A2GM,shapesGPA
    from labblouin.GM2PDB import write_coor, parsePDB,AtomInfo,WritePDB,Coord
    from glob import glob as G
    
    if ('*' in traj) or (type(traj) == list):
        # agglomerate all pdb files into a model file
        files = []
        [files.extend(G(t)) for t in traj]
        pdb = PDBnet.PDBstructure(ismodel=True)
        for f in range(len(files)):
            if not isdir('_input'): mkdir('_input')
            copy(files[f], '_input/')
            modelindex=f
            p = PDBnet.PDBstructure(files[f])
            if p.ismodel:
                for i in range(len(p.models.keys())):
                    key=p.models.keys()[i]
                    pdb.AddModel(i, p.models[key])
            else:
                for i in range(len(p.chains.keys())):
                    key=p.chains.keys()[i]
                    pdb.AddModel(i, p.chains[key]) 
    elif isfile(traj):
        #get pdbstructure from PDBnet
        pdb = PDBnet.PDBstructure(filein=traj)    
    #write the AA sequence for each model
    #create the fasta file includes AA sequences for all the models
    if not fasta:
        seq = pdb.ModelAsFASTA(1)
        if path.isfile('dummy.fasta'):remove('dummy.fasta')
        pdb.makefasta()
        '''fout=open("dummy.fasta",'w')
        for s in xrange(len(pdb.models)):
            fout.write('>model%d\n%s\n'%(s,seq))
        fout.close()'''
        fas = 'dummy.fasta'
    else: fas = fasta
    labels,coords,landmarkinfo = pdb.gm(fas,typeof='matrix',atomtype=atomtype)
    with open(prefix+'.landmarks','w') as F:
        for k, v in landmarkinfo.iteritems():
            F.write('>%s\n%s\n'%(k,''.join(v)))
    row, col = coords.shape
    k = col/dim
    arr = Rmatrix2Rarray(array2R(coords), dim=3)
    rot = shapesGPA(arr,scale=False)
    mat = A2GM(rot, dim=3)
    #transform every number from possible scientific notation to float notation
    if not fas: out = np.column_stack((range(len(labels)),mat))
    else: out = np.column_stack((labels,mat))
    fn = prefix+".gm"
    coord_df = pd.DataFrame(out)
    coord_df.to_csv(fn, sep=';', index = False, header=False)
    #np.savetxt(fn,out,delimiter=';',fmt='%s')
    if atomtype == 'all': 
        write_coor(parsePDB(pdb), prefix)
        lines = AtomInfo(prefix)
        d = Coord(fn)
        WritePDB(prefix, lines, d)


# Main Function

def main(arg, prefix='alignment',log=False,procrustes=None,alignr='matt',executable=None,nofasta=False,multi=0,
         MD=False,scores='GDT',tar=False,keep=False,split=False,reference=list(),cleanInput=False,cleanOutput=False,
         debug=False,scop=None,scopcache=None,atomtype='centroid',optimize=False,quick=0,gm=False):  
    
    # Setup logfile.
    if path.isfile('%s.log' % (prefix)):
        # Get current date/time as string for filename.
        ti = datetime.now().strftime('%Y_%m_%d_%H%M%S')
        newname = prefix + '.log.%s.old' % (ti)
        rename('%s.log' % (prefix),newname)
    log = logfile('%s.log' % (prefix),log)#options.log)

    # Greeter.
    log.write('ABeRMuSA: Approximate Best Reference Multiple ' + \
              'Structure Alignment ver. %s\n' % (VERSION),
              progress=False)        
    
    if procrustes:#options.procrustes:
        if path.isfile('%s.fasta'%(prefix)):fasta=prefix+'.fasta'
        else: fasta=None
        # Do the procrustes bit
        line = 'Executing partial Procrustes superimposition instead of traditional structural alignment\n'
        log.write(line,silent=False)
        Procrustes(prefix, arg,fasta=fasta,atomtype=procrustes,log=log) 
    else:  
        # Determine what binary/command to execute; see if recognized.
        cmd = alignr.lower()#options.aligner.lower()
        if cmd not in PLUGINS:
            # Given command/executable (aligner) is not supported.
            log.write('No plugin found under the name \'%s\'.' % (cmd),progress=False, status='ERROR')
            exit(127)
        aln = globals()[cmd]
    
        # Establish plugin rules, executable processing, multiprocessor support.
        if executable: exe = executable#if options.executable: exe = options.executable
        else: exe = aln.default_exe
        if not isCmd(exe):
            # Given command/executable not on system.
            log.write("'%s' is not present in your environment path as " % (exe) + \
                      'an executable. Please check your system path.',progress=False, status='ERROR')
            exit(126)
        elif not isCmd('muscle') and not nofasta:#options.nofasta:
            log.write("The executable 'muscle' was not found in your environment path.",
                progress=False, status='ERROR')
            exit(2)
        #exe = exewrapper(prefix,exe,aln,log,uniq=int(options.multi),ismodel=options.MD,ver=VERSION)
        exe = exewrapper(prefix,exe,aln,log,uniq=int(multi),ismodel=MD,ver=VERSION)
        if scores != None:#if options.scores != None:
            scores = scores.split(',') # Scores to do.
            for score in scores: exe.addScoreToDo(score.strip())
        log.write('Using \'%s\' as an aligner (plugin: \'%s\').' % (exe.cmd,cmd))
    
        # Load in files and folders.
        filelist, ref, t = [], None, T.getTime()    
        refwr = referenceCollection(reference) # Encapsulate references.
        acquireFiles(arg,filelist,log,refwr,cleanInput,split,MD)
    
        # Check files loaded for errors.
        if len(filelist) == 0:
            opts.print_help()
            exit(2)
        elif (len(filelist) == 1) and not MD:
            log.write('Only 1 argument <%s> was provided (need >= 2).' % (filelist[0]), status='ERROR')
            exit(2)
        if refwr.anyNotFound():
            if not MD:
                for ref in refwr.iterNotFound():
                    log.write('Reference <%s> was not in arguments. Adding.' % (ref), status='NOTE')
                    filelist.append(ref)   
    
        # Report success on loading files.
        log.write('%d protein structure files loaded (%ds).' % (len(filelist),T.getTime()-t))
    
        # Write an XML record of this run if debug mode is enabled.
        if debug:
            # Outline entire process as an "alignment".
            xml = XMLfile('%s.xml' % (prefix),'alignment') 
            # Record software's version number in XML file. 
            xml.root.set('version',VERSION)    
    
        # Setup SCOP and grab dissimilar PDBs if necessary.
        if scop and scopcache:
            if not path.isdir(PDB_CACHE): mkdir(PDB_CACHE)
            scinst = scop.scopHierarchy(scopcache)
            log.write('Acquiring SCOP metadata <cache %s>...' 
                      % (IO.getFolderName(scopcache)))
            scinst.populateHierarchy()
            _, pdbli = scinst.getDissimilar(scop)
            if pdbli:
                scoppdbs = []
                log.write('Downloading %d PDBs (RCSB)...' % (len(pdbli)+1))
                for pdb in pdbli:
                    log.writeTemporary('Downloading PDB <%s> to cache (%s)...' % (
                        pdb,PDB_CACHE))
                    orig = pfam.grabPDBFile(pdb,PDB_CACHE)
                    # Determine first chain.
                    p = PDBnet.PDBstructure(orig)
                    frstch = p.orderofchains[0]
                    log.writeTemporary('Extracting from <%s> chain %s...' % (pdb,frstch))
                    chai = pfam.extractPDBChain(orig,frstch,orig+'_%s' % (frstch))
                    scoppdbs.append(chai)
                exe.scpdbs = scoppdbs    
    
        # Handle quick reference search.
        qui = int(quick)
        if refwr.hasReference():
            quick = None
            if debug: xml.root.set('mode','guided')    
            log.write('Reference(s) were specified: %s.' % (', '.join(refwr.references)))
        elif quick:
            quick = quickref(filelist,qui,random=(qui >= 1))
            if debug:
                xml.root.set('mode','quick')
                xml.root.set('iterations','%d' % (qui))
        else:
            quick = None
            if debug: xml.root.set('mode','exhaustive')
            log.write('Performing an exhaustive search for best reference.')
    
        # Alert the user to what input was given and record to XML file.
        if debug:
            for f in filelist: xml.add(xml.root,'file',('xml',f),('folder',IO.getFileName(f)))
        log.write('Input (%d): %s' % (len(filelist),', '.join(filelist)),silent=True)
    
        # Set up estimated time remaining (timer).
        if ref: log.setTotalNum(2*len(filelist))
        elif quick: log.setTotalNum(2*quick.numiters*3*len(filelist))
        else: log.setTotalNum(2*((len(filelist))**2))
    
        # Perform alignment and parse output.
        log.write('-'*20 + 'Running Alignment Commands' + '-'*20)
        refdict, folders = aligner.run(filelist,log,ref=refwr,exe=exe,quick=quick)
        log.write('-'*20 + 'Parsing Alignment Output  ' + '-'*20)
        bestref, reffldr = aligner.parse(filelist,refdict,log)
    
        # See if parsing succeeded.
        if not bestref and not reffldr:
            log.writeElapsedTime()
            exit(1)
    
        # Write best reference to XML file.
        if debug:
            xml.add(xml.root,'reference',('xml',bestref[0]),('folder',reffldr),
                    ('rank','%.5f' % (bestref[1])))
    
        # Write alignment to PDB, FASTA, GM, and landmark files.
        log.write('Consolidating alignment into a set of single files...')
        m = multipleAlignment(
            filelist,prefix,bestref,reffldr,log,exe,atomtype,
            curate=cleanOutput,optimize=optimize,MD=MD,
            fasta=(not nofasta),gm=gm)
        status = m.construct()
    
        # See if GM file writing was successful.
        if not status:
            log.write('Process completed, but no PDB/FASTA/GM/landmark files were written.')
            log.writeElapsedTime()
            exit(1)
    
        # Tar all folders.
        if tar:
            tar = tarArchive(prefix)
            log.write('Compressing and archiving all reference folders (to %s)...' % (
                tar.fname))
            for fo in folders: tar.add(fo)
            tar.close()
            status.extend(tar.fname)     
    
        # Finalize and report success.
        log.updateTimer(log.totalnum)
        log.write('Process completed successfully.') 
        log.write('Final files included: %s.' % (', '.join(status)))
        if not debug and not keep:
            # Delete all reference folders if debug mode is not enabled and writing successful.
            log.write('Cleaning all reference trial folders (to keep these, use -k or enable debug mode)...', 
                status='NOTE')
            for fo in folders: rmDir(fo)
        else:
            log.write('Refraining from removing reference trial folders. Folders include: %s.' % 
                      (', '.join(folders)), status='NOTE')       
        log.writeElapsedTime()

###################################################################################################

# Set up option parsing.

opts = optparse.OptionParser(usage='%prog [options] file1/folder1 [file2/folder2]...')
opts.add_option('--log', '-l', action='store_true',default=False,
                help='Whether or not to write a logfile. Default: False.')
opts.add_option('--aligner', '-a', default='matt',dest='alignr',
                help='Specify the pairwise aligner to use. Only those aligners that are supported (possess a plugin) will be accepted. Supported: ' + ', '.join(PLUGINS) + '. Default: matt.')
opts.add_option('--reference','-r',action='append',type='string',default=list(),
                help='Specify a particular reference PDB file (by path) if necessary. Multiple references can be specified by prefixing this option before them. Default: no references are specified.')
opts.add_option('--quick','-q', default=0,
                help='If a reference is not provided and this number is greater than 0, use the quick reference search method instead of a full exhaustive search and iterate that number of times. Default: 0 (disabled; full exhaustive search).')
opts.add_option('--atomtype','-c',action='store_true',default='centroid',
                help='Whether or not to consider alpha carbons when writing to GM. Default: centroid.')
opts.add_option('--executable', '-e', default=None,
                help='If the intended executable for the aligner given is not on path, it can be specified here. If left blank, will use default value, dependent on plugin.')
opts.add_option('--prefix','-p', default='alignment',
                help='If a specific prefix is necessary, it can be provided here. Will default to alignment.')
opts.add_option('--optimize','-o',action='store_true',default=False,
                help='Whether or not to optimize during the writing of the GM file (by number of  landmarks). Feature still in testing. Default: False.')
opts.add_option('--scores','-s',default='GDT',
                help='Specify a scoring method you would like to be used on alignment. Default: GDT. Scoring methods include: %s. Multiple scores can be specified ' % (
                    ', '.join(SCORE_TYPES)) + 'by separating by commas.')
opts.add_option('--scop','-y',default=None,
                help='Specify a SCOP superfamily in the case of a small residue if encountered. Ensures accuracy in the case of statistical failure of RRMSD measure.'+\
                ' Requires scopcache to also be specified. Default: None.')
opts.add_option('--scopcache','-z',default=None,
                help='Specify the location of a SCOP cache locally if necessary. Default: None.')
opts.add_option('--multi','-m', default=0,
                help='Whether or not to perform execution on a multiprocessor platform, e.g. Fester. Default: 0. Anything greater than 0 will imply the use of a Grid Engine and will specify the number of cores necessary.')
opts.add_option('--nofasta','--nf',action='store_true',default=False,
                help='When enabled, will not consolidate all alignments into a single FASTA file, nor will it write either the GM or landmark file since they require a FASTA file to be written.')
opts.add_option('--cleanInput','--cinput', action='store_true', default=False,
                help='Whether or not to clean/curate input PDB files. Default: False.')
opts.add_option('--cleanOutput','--coutput', action='store_true', default=False,
                help='Whether or not to clean/curate output PDB files. Default: False.')
opts.add_option('--split','-x', action='store_true', default=False,
                help='Whether or not to split separate chains into separate files. Default: False.')
opts.add_option('--tar','-t', action='store_true', default=False,
                help='Whether or not to compress all output folders for reference PDBs into a compressed tarfile named by prefix. Can take a long time if set size is large.')
opts.add_option('--MD','-M', action='store_true', default=False,
                help='Whether or not to assume that the single file passed is a trajectory. Effectively this will split the trajectory, align every model.')
opts.add_option('--keep','-k',action='store_true',default=False,
                help='Retain reference folders after completion of the alignment and do not delete them.')
opts.add_option('--debug','-D',action='store_true',default=False,
                help='Enable this to output debug information and metadata about executed alignments into XML file(s). Also retains reference folders after completion.')
opts.add_option('--procrustes','-P',action='store',default=None,
                help='Perform partial procrustes superimposition instead of the regular alignment. Only suitable for md simulations for now. Values should be centroid, CA or all')
opts.add_option('--Pairwise2GM','--gm',action='store_true',default=False, dest='gm',
                help='If used in addition of nofasta (--nf) a GM file will be written based on the homologous residues in all pairwise alignments.')

# If not imported.

if __name__ == "__main__":
    # Option Parsing.
    options, arg = opts.parse_args()      
    main(arg,**vars(options))
