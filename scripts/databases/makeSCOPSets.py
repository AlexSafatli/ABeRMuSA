#!/usr/bin/env python

# makeSCOPSets.py
# -------------------------
# Jul 17, 2013; Alex Safatli
# -------------------------
# Walk through SCOP Database
# and create PDB sets.

import sys, os, shutil, random
from utils import scop, IO, pfam
from utils.logfile import XMLfile

if len(sys.argv) != 7: 
    print 'usage: %s SCOP_METADATA_PATH PDB_CACHE NUM_SETS MIN_NUM MAX_NUM STEP_NUM' % (sys.argv[0])
    exit()

scpath = sys.argv[1]
print '----> SCOP Metadata Path: %s' % (scpath)
pdcache = sys.argv[2]
print '----> PDB File Cache Path: %s' % (pdcache)
numsets = int(sys.argv[3])
print '----> Number of Sets to Generate per Range: %d' % (numsets)
maxnum = int(sys.argv[4])
minnum = int(sys.argv[5])
stenum = int(sys.argv[6])
numrange = range(maxnum,minnum,stenum)
print '----> Range for Number of Items per Set to Generate: %s' % (numrange)
usedsuperfamilies = []

# Get SCOP metadata.
scinst = scop.scopHierarchy(scpath)
print '====> Acquiring SCOP metadata <cache %s>' % (IO.getFolderName(scpath))
scinst.populateHierarchy()
print '====> Acquiring PDB files and set data'

# For NUM_SETS, begin to generate a set with NUM_IN_SET items that are
# similar and NUM_IN_SET items that are dissimilar.
for i in range(0,numsets):
    
    for numper in numrange:
    
        # Get lists ready for PDB paths.
        dis_paths = []
        sim_paths = []
        
        # Get a (unique) superfamily.
        print '----> Set number: %d' % (i)
        sfs = scinst.superfamilies.keys()
        sf = None
        while (sf == None or sf in usedsuperfamilies):
            rk = random.randint(0,len(sfs))
            sf = sfs[rk]
        usedsuperfamilies.append(sf)
        print 'Superfamily: %s' % (sf)
        print 'Number in set: %d' % (numper)
    
        # Get NUM_IN_SET dissimilar PDBs.
        _, pdbli_d = scinst.getDissimilar(sf,numper)
    
        # Get NUM_IN_SET similar PDBs.
        _, pdbli_s = scinst.getSimilar(sf,numper)
        
        # Download necessary PDBs, first checking cache.
        for pdb in pdbli_d:
            print 'Downloading %s <dissimilar>...' % (pdb)
            o = pfam.grabPDBFile(pdb,pdcache)
            dis_paths.append(o)
        for pdb in pdbli_s:
            print 'Downloading %s <similar>...' % (pdb)
            o = pfam.grabPDBFile(pdb,pdcache)
            sim_paths.append(o)
            
        # Save set to file system.
        if os.path.isfile('%s-%d.xml' % (sf,numper)):
            repl = raw_input('Replace %s.xml with new superfamily data (Y/N)? ')
            ans  = (repl.upper() == 'Y')
            if not ans:
                print 'Superfamily data scrapped. No file was replaced.'
                continue
        xml = XMLfile('%s-%d.xml' % (sf,numper),'superfamily')
        for pdb in dis_paths: xml.add(xml.root,'dissimilar',('xml',pdb))
        for pdb in sim_paths: xml.add(xml.root,'similar',('xml',pdb))
    