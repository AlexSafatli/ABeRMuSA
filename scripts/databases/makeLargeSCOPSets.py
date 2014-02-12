#!/usr/bin/env python

# makeLargeSCOPSets.py
# -------------------------
# Aug 13, 2013; Alex Safatli
# -------------------------
# Prepares a collection of
# SCOP set files for use in
# traversing sections of the
# SCOP database.
#
# Walk through SCOP Database
# and create 1 PDB set ranging
# from a minimum to a maximum
# by an exponential function;
# uses only similar superfamily
# PDBs for any given set and extracts
# only the first chain from any PDB.

import sys, os, shutil, random
from utils import scop, IO, pfam, PDBnet
from utils.logfile import XMLfile

if len(sys.argv) != 6: 
    print 'usage: %s SCOP_METADATA_PATH PDB_CACHE MIN_NUM MAX_NUM EXP_BASE' % (sys.argv[0])
    exit()

scpath = sys.argv[1]
print '----> SCOP Metadata Path: %s' % (scpath)
pdcache = sys.argv[2]
print '----> PDB File Cache Path: %s' % (pdcache)
minnum = int(sys.argv[3])
maxnum = int(sys.argv[4])
basenm = int(sys.argv[5])
numrange = []

# Generate range.
num = 0
exp = 0
while (num <= maxnum):
    if (num != 0): numrange.append(num)
    num = basenm**exp+minnum
    exp += 1
    
print '----> Range for Set to Generate: %s' % (numrange)
usedsuperfamilies = []

# Get SCOP metadata.
scinst = scop.scopHierarchy(scpath)
print '====> Acquiring SCOP metadata <cache %s>' % (IO.getFolderName(scpath))
scinst.populateHierarchy()
print '====> Acquiring PDB files and set data'

# Get a (unique) superfamily.
sfs = scinst.superfamilies.keys()
sf = None
while (sf == None or sf in usedsuperfamilies 
       or len(scinst.superfamilies[sf].getPDBs()) < maxnum):
    rk = random.randint(0,len(sfs))
    sf = sfs[rk]
usedsuperfamilies.append(sf)
print '====> Superfamily: %s' % (sf)
_, pdbli = scinst.getSimilar(sf,numrange[-1])

for numper in numrange:
    
    # See if already formed.
    if os.path.isfile('%s-%d.xml' % (sf,numper)):
        repl = raw_input('Replace %s.xml with new superfamily data (Y/N)? ')
        ans  = (repl.upper() == 'Y')
        if not ans:
            print 'No file was replaced.'
            continue  
    
    # Get lists ready for PDB paths.
    sim_paths = []
    print '----> Number in set: %d' % (numper)
    pdbli_s = pdbli[:numper]

    # Download necessary PDBs, first checking cache.
    for pdb in pdbli_s:
        on = os.path.join(pdcache,pdb+'.pdb')
        if not os.path.isfile(on):
            print 'Downloading %s <similar>...' % (pdb)
            o = pfam.grabPDBFile(pdb,pdcache)
        else: o = on
        pdinst = PDBnet.PDBstructure(o)
        frstch = pdinst.orderofchains[0]
        fn = os.path.join(
            pdcache,IO.getFileName(o)+'_%s.pdb' % (frstch))
        if not os.path.isfile(fn):
            print 'Extracting %s from %s...' % (frstch,pdb)
            p = pfam.extractPDBChain(o,frstch,fn)
        else: p = fn
        sim_paths.append(p)
    
    # Save to disk.
    xml = XMLfile('%s-%d.xml' % (sf,numper),'superfamily')
    for pdb in sim_paths: xml.add(xml.root,'similar',('xml',pdb))
    