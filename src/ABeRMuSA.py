#!/usr/bin/env python

''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca

+----------------------------------------------+
| ABeRMuSA.py < Automatic Pairwise Alignment > |
+----------------------------------------------+
May 6th, 2013; Alex Safatli; Complete refactoring of Kyle Nguyen code (2012).

Automatic pairwise alignment capable of being extended to any pairwise alignment executable. By default, uses MATT. ABeRMuSA stands for "Approximate Best Reference Multiple Structure Alignment" method.

Input:   PDB files, folders of PDB files.
Options: See help menu (--help, -h). '''

# Imports

import optparse, glob, sys, aligner, tarfile
from utils import homology, pfam, scop, PDBnet, IO
from utils.logfile import logfile, XMLfile, timer as T
from os import path, mkdir, rename
from datetime import datetime
from exewrapper import exewrapper, exeExists as isCmd
from quickref import quickref
from writer import multipleAlignment
from scoring import SCORE_TYPES
from plugins import * # get all plugins

# Constants, Initialization

SCRIPT_FOLDER = path.split(path.realpath(__file__))[0]
PLUGIN_FOLDER = path.join(SCRIPT_FOLDER,'plugins')
PDB_CACHE     = path.join(SCRIPT_FOLDER,'pdbcache')
PLUGIN_PYS    = glob.glob(path.join(PLUGIN_FOLDER,'*.py'))
PLUGINS       = [path.split(x)[-1].strip('.py') for x in PLUGIN_PYS \
                 if not x.endswith('__init__.py')]
VERSION       = '0.5.0'
PDB_ALLOW     = ['pdb','ent','atm']
PDB_FOLDER    = '_input'

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
    if path.isfile(fn): return fn
    towrite = []
    o = open(fi)
    l = o.readlines()
    o.close()
    for li in l:
        if li.startswith('ATOM'): towrite.append(li)
    o = open(fn,'w')
    o.write(''.join(towrite))
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
    

def handleFile(fi,log,refw,clean=False,split=False):
    
    ''' Handle a file during startup; do cleaning, etc. '''
    
    # Is it the reference?
    isref = False
    if refw.ref == fi:
        refw.found = True
        isref = True
    
    # Need only ATOM lines.
    fn = IO.getFileName(fi) # Get filename.
    af = path.join(PDB_FOLDER,fn+'.pdb')
    fi = stripAtomsFromPDB(fi,af)
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
                log.write('Extracted <%s> (%s) to <%s>.' % (fn,ch,cf))
                multi.append(cf)
            return multi
    
    # Return the file.
    return [fi]
        

def acquireFiles(arg,fl,log,refw,clean=False,split=False):
    
    ''' Populate an input filelist with alignable items from 
    the arguments. '''
    
    # Enumerate over arguments.
    for fi in arg:
        # Is a directory?
        if path.isdir(fi):
            # Enumerate allowable items in directory.
            darg = enumerateFilesInFolder(fi)
            acquireFiles(darg,fl,log,refw,clean,split)
        # Is a file?
        elif path.isfile(fi):
            fl.extend(handleFile(fi,log,refw,clean,split))
        # Is not file BUT is reference?
        elif fi == refw.ref:
            log.write('ERROR; <%s> assigned as reference; does not exist.' % (
                fi))
            exit(1)
        # Not file or folder.
        else: log.write('WARNING; <%s> not a file or folder.' % (fi))

# Main Function

def main(options,arg):  
    
    # Determine a prefix for operation.
    prefix = options.prefix
    if not prefix: prefix = 'alignment'     
    
    # Setup logfile.
    log = logfile('%s.log' % (prefix),options.log)    
    
    # Greeter.
    log.write('ABeRMuSA: Approximate Best Reference Multiple ' + \
              'Structure Alignment ver. %s\n' % (VERSION))      
    
    # Determine what binary/command to execute; see if recognized.
    t = T.getTime()
    cmd = options.aligner.lower()
    if cmd not in PLUGINS:
        # Given command/executable (aligner) is not supported.
        log.write('%s is not a pairwise aligner currently supported ' + \
                  'by this software. No plugin found.' % (cmd))
        exit(2)
    aln = globals()[cmd]     
    
    # Establish plugin rules, executable processing, multiprocessor support.
    if options.executable: exe = options.executable
    else: exe = aln.default_exe
    if not isCmd(exe):
        # Given command/executable not on system.
        log.write('%s is not present in your environment ' + \
                  'path as an executable. Please check your system path.' % (exe))
        exit(2)    
    exe = exewrapper(prefix,exe,aln,log,uniq=int(options.multi))
    if options.scores != None:
        scores = options.scores.split(',') # Scores to do.
        for score in scores: exe.addScoreToDo(score.strip())
    log.write('Plugin %s loaded (%ds).' % (cmd,T.getTime()-t))    
    
    # See if folder for PDBs has been made; otherwise, make it.
    if not path.isdir(PDB_FOLDER): mkdir(PDB_FOLDER)
    log.write('Created folder %s to store PDB files for input to alignment.' % (
        PDB_FOLDER))
    
    # Load in files and folders.
    t = T.getTime()
    refwr = refWrapper(options.reference)
    filelist, ref = [], None
    acquireFiles(arg,filelist,log,refwr,
                 options.clean,options.split)
    
    # Check files loaded for errors.
    if len(filelist) == 0:
        opts.print_help()
        exit(2)
    elif len(filelist) == 1:
        log.write('ERROR; Only 1 argument <%s> was provided (need >= 2).' 
                  % (filelist[0]))
        exit(2)
    if not options.reference: pass
    elif not refwr.found:
        log.write('WARNING; Reference <%s> was not in arguments. Adding.' 
                  % (options.reference))
        filelist.append(options.reference)
        ref = options.reference
    else: ref = refwr.ref    

    # Report success on loading files.
    log.write('All (%d) files loaded (%ds).' % (len(filelist),T.getTime()-t))
    
    # Write an XML record of this run.
    if path.isfile('%s.xml' % (prefix)):
        ti = datetime.now().strftime('%Y_%m_%d_%H%M%S') # Get time as string for filename.
        log.write('WARNING; An XML file already exists with prefix name <%s>.' % (prefix))
        log.write('Moving previous XML file to %s.%s.old' % ('%s.xml' % (prefix),ti))
        rename('%s.xml' % (prefix),prefix + '.xml.%s.old' % (ti))
    xml = XMLfile('%s.xml' % (prefix),'alignment') # Outline entire process as an "alignment".
    xml.root.set('version',VERSION) # Record software's version number in XML file.    
    
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
            log.writeTemporary('Downloading %d PDBs (RCSB)...\n' % (
                len(pdbli)+1))
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
        xml.root.set('mode','guided')    
    elif options.quick:
        quick = quickref(filelist,qui,random=(qui >= 1))
        xml.root.set('mode','quick')
        xml.root.set('iterations','%d' % (qui))
    else:
        quick = None
        xml.root.set('mode','exhaustive')
    
    # Alert the user to what input was given and record to XML file.
    for f in filelist: xml.add(xml.root,'file',('xml',f),('folder',IO.getFileName(f)))
    for f in filelist: log.write('Input: %s' % (f))    
    
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
    xml.add(xml.root,'reference',('xml',bestref[0]),('folder',reffldr),
            ('rank','%.5f' % (bestref[1])))
    
    # Write alignment to PDB, FASTA, GM, and landmark files.
    log.write('Consolidating alignment...')
    m = multipleAlignment(
        filelist,prefix,bestref,reffldr,log,exe,options.alphaC,
        optimize=options.optimize)
    status = m.construct()
    
    # See if GM file writing was successful.
    if not status:
        log.write('Process completed; no PDB/FASTA/GM/landmark files written.')
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
                help='Specify the pairwise aligner to use. Only ' +\
                'those aligners that are supported (possess a plugin) will be accepted. Supported: '+\
                ', '.join(PLUGINS) + '. Default: matt.')
opts.add_option('--reference','-r',default=None,
                help='Specify a particular reference PDB file if necessary. Default: None.')
opts.add_option('--quick','-q', default=0,
                help='If a reference is not provided and this number is greater than 0, '+\
                'use the quick reference search method instead of a full exhaustive search and '+
                'iterate that number of times. Default: 0 (full exhaustive search).')
opts.add_option('--alphaC','-c',action='store_true',default=False,
                help='Whether or not to consider alpha carbons when writing to GM. Default: False.')
opts.add_option('--executable', '-e', default=None,
                help='If the intended executable for the aligner given is not on '+\
                'path, it can be specified here. If left blank, will use default value, dependent on plugin.')
opts.add_option('--prefix','-p', default=None,
                help='If a specific prefix is necessary, it can be provided here. Will default'+\
                ' to alignment.')
opts.add_option('--optimize','-o',action='store_true',default=False,
                help='Whether or not to optimize during the writing of the GM file (by number of'+\
                ' landmarks). Feature still in testing. Default: False.')
opts.add_option('--scores','-s',default='TMscore',
                help='Specify a scoring method you would like to be used on alignment. Default:' +\
                ' TMscore. Scoring methods include: %s. Multiple scores can be specified ' % (
                    ', '.join(SCORE_TYPES)) + 'by separating by commas.')
opts.add_option('--scop','-y',default=None,
                help='Specify a SCOP superfamily in the case of a small residue if '+\
                'encountered. Ensures accuracy in the case of statistical failure of RRMSD measure.'+\
                ' Requires scopcache to also be specified. Default: None.')
opts.add_option('--scopcache','-z',default=None,
                help='Specify the location of a SCOP cache locally if necessary. Default: None.')
opts.add_option('--multi','-m', default=0,
                help='Whether or not to perform execution on a multiprocessor platform, e.g. Fester. '+\
                'Default: 0. Anything greater than 0 will imply the use of a Grid Engine and will specify '+\
                'the number of cores necessary.')
opts.add_option('--clean','-n', action='store_true', default=False,
                help='Whether or not to clean/curate input PDB files. Default: False.')
opts.add_option('--split','-x', action='store_true', default=False,
                help='Whether or not to split separate chains into separate files. Default: False.')
opts.add_option('--tar','-t', action='store_true', default=False,
                help='Whether or not to compress all output folders for reference PDBs into a compressed'+\
                ' tarfile named by prefex. Can take a long time if set size is large.')

# If not imported.

if __name__ == "__main__":
    # Option Parsing.
    options, arg = opts.parse_args()      
    main(options, arg)
