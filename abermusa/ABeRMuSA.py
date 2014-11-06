#!/usr/bin/env python

''' 
ABeRMuSA (Automatic Pairwise Alignment); complete refactoring of Kyle Nguyen code (2012). 
Automatic pairwise alignment capable of being extended to any pairwise alignment executable. 
By default, uses MATT. ABeRMuSA stands for "Approximate Best Reference Multiple Structure Alignment"
method.

Input:   PDB files, folders of PDB files.
Options: See help menu (--help, -h). 
'''

# Date:   May 6 2013
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca
# Contributions made by Jose Sergio Hleap (2014)

# Imports

import optparse, glob, sys, aligner, tarfile
from labblouin import homology, pfam, scop, PDBnet, IO
from labblouin.logfile import logfile, XMLfile, timer as T
from os import path, mkdir, rename
from shutil import rmtree as rmDir
from datetime import datetime
from exewrapper import exewrapper, exeExists as isCmd
from quickref import quickref
from writer import multipleAlignment
from scoring import SCORE_TYPES
from plugins import * # get all plugins

# Constants, Initialization

SCRIPT_FOLDER = path.split(path.realpath(__file__))[0]
PLUGIN_FOLDER = path.join(SCRIPT_FOLDER,'plugins')
PLUGIN_PYS    = glob.glob(path.join(PLUGIN_FOLDER,'*.py'))
PLUGINS       = [path.split(x)[-1].strip('.py') for x in PLUGIN_PYS \
                 if not x.endswith('__init__.py')]
VERSION       = '0.5.1.3'
PDB_ALLOW     = ['pdb','ent','atm']
PDB_FOLDER    = '_input'
PDB_CACHE     = '_scop'

__author__    = 'Alex Safatli'
__version__   = VERSION

# Auxiliary Functions, Classes

class tarWrapper:
    def __init__(self,prefix):
        self.prefix = prefix
        self.fname  = prefix + '.tar.gz'
        self.tar    = tarfile.open(self.fname,"w:gz")
    def add(self,fi):
        self.tar.add(fi)
        return fi
    def close(self):
        self.tar.close()

class refWrapper:
    def __init__(self,ref):
        self.ref = ref
        self.found = False

def stripAtomsFromPDB(fi,fn):
    
    ''' Get atomlines from fi file and write to fn. '''
    
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

    # If there is a reference pickle file in path already.
    refpi = path.join(path,'ref.pickl')
    if path.isfile(refpi): return []

    # Go over all allowable extensions.
    for ext in PDB_ALLOW:
        lookf = path.join(path,'*.' + ext)
        enum  = glob.glob(lookf)
        files.extend(enum)

    return files

def splitTraj(filename,log):
    
    ''' Split a trajectory PDB file (from a molecular dynamics simulation) into individual files. '''
    
    multi = []
    with open(filename) as F:
        model=''
        for line in F:
            if line.startswith('MODEL'):
                bl = line.strip().split()
                num = int(bl[1])
                mod = num-1
                if mod == 0: continue
                else:
                    cf = path.join(PDB_FOLDER,'Model%d.pdb'%(mod))
                    if not path.isfile(cf):
                        log.writeTemporary('Extracted <%s> (Model %s) to <%s>.' % (filename,mod,cf))
                        fout=open(cf,'w')
                        fout.write(model)
                        fout.close()
                        model=''	
                        multi.append(cf)
                    else:
                        log.writeTemporary('File Model%s.pdb already in input folder.' % (mod))
                        model=''	
                        multi.append(cf)                        
            elif line.startswith('ATOM'): model+=line
            else: continue    
    return multi

def handleFile(fi,log,refw,clean=False,split=False,MD=False):

    ''' Handle a file during startup; do cleaning, etc. '''

    # Is it the reference?
    isref = False
    if refw.ref == fi:
        refw.found = True
        isref = True

    # Need only ATOM lines.
    fn = IO.getFileName(fi) # Get filename.
    if not MD:
        af = path.join(PDB_FOLDER,fn+'.pdb')
        fi = stripAtomsFromPDB(fi,af)
    else: fi = fn
    if isref: refw.ref = fi

    # Translate to PDBnet structure if necessary.
    if (clean or split):
        try: p = PDBnet.PDBstructure(fi)
        except:
            log.write('WARNING; <%s> could not be parse as PDB.' % (fn))
            if isref:
                log.write('ERROR; <%s> is provided reference.' % (fn))
                exit(1)
            return []

    # See if needs to be cleaned.
    if (clean and not p.CheckComplete()):
        cf = path.join(PDB_FOLDER,fn+'-c.pdb')
        if not path.isfile(cf): homology.completePDB(fi,cf)
        fi = cf
        if isref: refw.ref = fi    

    # Multiple chains?
    if split:
        chains = p.orderofchains
        numch  = len(chains)
        if (numch > 1):
            # Split it.
            log.writeTemporary('NOTE; <%s> had multiple chains (%d). Extracting.' % (
                fn,numch))
            multi = []
            for ch in chains:
                cf = path.join(PDB_FOLDER,fn+'_%s.pdb' %(ch))
                if not path.isfile(cf): pfam.extractPDBChain(fi,ch,cf)
                log.writeTemporary('Extracted <%s> (%s) to <%s>.' % (fn,ch,cf))
                multi.append(cf)
            return multi

    # Multiple models?
    if MD:
        # Split it.
        log.writeTemporary('NOTE; <%s> had multiple models. Extracting.' % (
            fn))
        multi = splitTraj(fn+'.pdb',log)
        refw.found = True
        refw.ref = path.join(PDB_FOLDER,refw.ref)
        return multi		

    # Return the file.
    return [fi]

def acquireFiles(arg,fl,log,refw,clean=False,split=False,MD=False):

    ''' Populate an input filelist with alignable items from 
    the arguments. '''

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
        elif fi == refw.ref:
            writeTemporary('ERROR; <%s> assigned was reference but does not exist.' % (
                fi))
            exit(1)
        # Not file or folder.
        else: lwriteTemporary('WARNING; <%s> not a file or folder.' % (fi))

# Main Function

def main(options,arg):  

    # Determine a prefix for operation.
    prefix = options.prefix
    if not prefix: prefix = 'alignment'     

    # Setup logfile.
    if path.isfile('%s.log' % (prefix)):
        ti = datetime.now().strftime('%Y_%m_%d_%H%M%S') # Get time as string for filename.
        rename('%s.log' % (prefix),prefix + '.log.%s.old' % (ti))    
    log = logfile('%s.log' % (prefix),options.log)    

    # Greeter.
    log.write('ABeRMuSA: Approximate Best Reference Multiple ' + \
              'Structure Alignment ver. %s\n' % (VERSION))      

    # Determine what binary/command to execute; see if recognized.
    cmd = options.aligner.lower()
    if cmd not in PLUGINS:
        # Given command/executable (aligner) is not supported.
        log.write('%s is not a pairwise aligner currently supported ' % (cmd)) + \
            'by this software. No plugin found.' 
        exit(2)
    aln = globals()[cmd]     

    # Establish plugin rules, executable processing, multiprocessor support.
    if options.executable: exe = options.executable
    else: exe = aln.default_exe
    if not isCmd(exe):
        # Given command/executable not on system.
        log.write('%s is not present in your environment path as ' % (exe) + \
                  'an executable. Please check your system path.')
        exit(2)    
    exe = exewrapper(prefix,exe,aln,log,uniq=int(options.multi),ismodel=options.MD,ver=VERSION)
    if options.scores != None:
        scores = options.scores.split(',') # Scores to do.
        for score in scores: exe.addScoreToDo(score.strip())
    log.write('Plugin %s loaded.' % (cmd))

    # See if folder for PDBs has been made; otherwise, make it.
    if not path.isdir(PDB_FOLDER): mkdir(PDB_FOLDER)
    log.write('Created folder %s to store PDB files for input to alignment.' % (
        PDB_FOLDER))

    # Load in files and folders.
    filelist, ref, t = [], None, T.getTime()    
    refwr = refWrapper(options.reference) # Encapsulate reference.
    acquireFiles(arg,filelist,log,refwr,options.cleanInput,options.split,options.MD)

    # Check files loaded for errors.
    if len(filelist) == 0:
        opts.print_help()
        exit(2)
    elif (len(filelist) == 1) and not options.MD:
        log.write('ERROR; Only 1 argument <%s> was provided (need >= 2).' 
                  % (filelist[0]))
        exit(2)
    if not options.reference: pass
    elif not refwr.found:
        if not options.MD:
            log.write('WARNING; Reference <%s> was not in arguments. Adding.' 
                      % (options.reference))
            filelist.append(options.reference)
            ref = options.reference
        else:
            ref = path.join(PDB_FOLDER,options.reference)
            if ref in filelist:
                refwr.found=True
                refwr.ref = ref
                
    else: ref = refwr.ref    

    # Report success on loading files.
    log.write('All (%d) files loaded (%ds).' % (len(filelist),T.getTime()-t))

    # Write an XML record of this run if debug mode is enabled.
    if options.debug:
        # Outline entire process as an "alignment".
        xml = XMLfile('%s.xml' % (prefix),'alignment') 
        # Record software's version number in XML file. 
        xml.root.set('version',VERSION)    

    # Setup SCOP and grab dissimilar PDBs if necessary.
    if options.scop and options.scopcache:
        if not path.isdir(PDB_CACHE): mkdir(PDB_CACHE)
        scinst = scop.scopHierarchy(options.scopcache)
        log.write('Acquiring SCOP metadata <cache %s>...' 
                  % (IO.getFolderName(options.scopcache)))
        scinst.populateHierarchy()
        _, pdbli = scinst.getDissimilar(options.scop)
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
    qui = int(options.quick)
    if ref:
        quick = None
        if options.debug: xml.root.set('mode','guided')    
    elif options.quick:
        quick = quickref(filelist,qui,random=(qui >= 1))
        if options.debug:
            xml.root.set('mode','quick')
            xml.root.set('iterations','%d' % (qui))
    else:
        quick = None
        if options.debug: xml.root.set('mode','exhaustive')

    # Alert the user to what input was given and record to XML file.
    if options.debug:
        for f in filelist: xml.add(xml.root,'file',('xml',f),('folder',IO.getFileName(f)))
    log.write('Input (%d): %s' % (len(filelist),', '.join(filelist)))

    # Set up estimated time remaining (timer).
    if ref: log.setTotalNum(2*len(filelist))
    elif options.quick: log.setTotalNum(2*quick.numiters*3*len(filelist))
    else: log.setTotalNum(2*((len(filelist))**2))

    # Perform alignment and parse output.
    log.write('-'*20 + 'Running Alignment Commands' + '-'*20)
    refdict, folders = aligner.run(filelist,log,ref=ref,exe=exe,quick=quick)
    log.write('-'*20 + 'Parsing Alignment Output  ' + '-'*20)
    bestref, reffldr = aligner.parse(filelist,refdict,log)

    # See if parsing succeeded.
    if not bestref and not reffldr:
        log.writeElapsedTime()
        exit(1)

    # Write best reference to XML file.
    if options.debug:
        xml.add(xml.root,'reference',('xml',bestref[0]),('folder',reffldr),
                ('rank','%.5f' % (bestref[1])))

    # Write alignment to PDB, FASTA, GM, and landmark files.
    log.write('Consolidating alignment into a set of single files...')
    m = multipleAlignment(
        filelist,prefix,bestref,reffldr,log,exe,options.alphaC,
        curate=options.cleanOutput,optimize=options.optimize,MD=options.MD)
    status = m.construct()

    # See if GM file writing was successful.
    if not status:
        log.write('Process completed but no PDB/FASTA/GM/landmark files were written.')
        log.writeElapsedTime()
        exit(1)

    # Tar all folders.
    if options.tar:
        tar = tarWrapper(prefix)
        log.write('Compressing and archiving all reference folders (to %s)...' % (
            tar.fname))
        for fo in folders: tar.add(fo)
        tar.close()
        status.extend(tar.fname)
    if not options.debug:
        # Delete all reference folders if debug mode is not enabled and writing successful.
        log.write('Removing all reference folders...')
        for fo in folders: rmDir(fo)
        
    # Finalize and report success.
    log.updateTimer(log.totalnum)
    log.write('Process completed successfully.')
    log.write('Associated folders output for reference trials include: %s.' % 
              (', '.join(folders)))    

    log.write('Final files include: %s.' % (', '.join(status)))
    log.writeElapsedTime()

###################################################################################################

# Set up option parsing.

opts = optparse.OptionParser(usage='%prog [options] file1/folder1 [file2/folder2]...')
opts.add_option('--log', '-l', action='store_true',default=False,
                help='Whether or not to write a logfile. Default: False.')
opts.add_option('--aligner', '-a', default='matt',
                help='Specify the pairwise aligner to use. Only those aligners that are supported (possess a plugin) will be accepted. Supported: ' + ', '.join(PLUGINS) + '. Default: matt.')
opts.add_option('--reference','-r',default=None,
                help='Specify a particular reference PDB file if necessary. Default: None.')
opts.add_option('--quick','-q', default=0,
                help='If a reference is not provided and this number is greater than 0, use the quick reference search method instead of a full exhaustive search and iterate that number of times. Default: 0 (disabled; full exhaustive search).')
opts.add_option('--alphaC','-c',action='store_true',default=False,
                help='Whether or not to consider alpha carbons when writing to GM. Default: False.')
opts.add_option('--executable', '-e', default=None,
                help='If the intended executable for the aligner given is not on path, it can be specified here. If left blank, will use default value, dependent on plugin.')
opts.add_option('--prefix','-p', default=None,
                help='If a specific prefix is necessary, it can be provided here. Will default to alignment.')
opts.add_option('--optimize','-o',action='store_true',default=False,
                help='Whether or not to optimize during the writing of the GM file (by number of  landmarks). Feature still in testing. Default: False.')
opts.add_option('--scores','-s',default='TMscore',
                help='Specify a scoring method you would like to be used on alignment. Default: TMscore. Scoring methods include: %s. Multiple scores can be specified ' % (
                    ', '.join(SCORE_TYPES)) + 'by separating by commas.')
opts.add_option('--scop','-y',default=None,
                help='Specify a SCOP superfamily in the case of a small residue if encountered. Ensures accuracy in the case of statistical failure of RRMSD measure.'+\
                ' Requires scopcache to also be specified. Default: None.')
opts.add_option('--scopcache','-z',default=None,
                help='Specify the location of a SCOP cache locally if necessary. Default: None.')
opts.add_option('--multi','-m', default=0,
                help='Whether or not to perform execution on a multiprocessor platform, e.g. Fester. Default: 0. Anything greater than 0 will imply the use of a Grid Engine and will specify the number of cores necessary.')
opts.add_option('--cleanInput','--cinput', action='store_true', default=False,
                help='Whether or not to clean/curate input PDB files. Default: False.')
opts.add_option('--cleanOutput','--coutput', action='store_true', default=False,
                help='Whether or not to clean/curate output PDB files. Default: False.')
opts.add_option('--split','-x', action='store_true', default=False,
                help='Whether or not to split separate chains into separate files. Default: False.')
opts.add_option('--tar','-t', action='store_true', default=False,
                help='Whether or not to compress all output folders for reference PDBs into a compressed tarfile named by prefex. Can take a long time if set size is large.')
opts.add_option('--MD','-M', action='store_true', default=False,
                help='Whether or not to assume that the single file passed is a trajectory. Effectively this will split the trajectory, align every model.')
opts.add_option('--debug','-D',action='store_true',default=False,
                help='Enable this to output debug information and metadata about executed alignments including XML files. Also retains reference folders.')

# If not imported.

if __name__ == "__main__":
    # Option Parsing.
    options, arg = opts.parse_args()      
    main(options, arg)
