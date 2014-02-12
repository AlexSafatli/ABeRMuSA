#!/usr/bin/env python

# makeSCOPSetForSuperfamily.py
# -------------------------
# Jul 17, 2013; Alex Safatli
# -------------------------
# Grab all families associated
# with a superfamily and all
# corresponding PDBs. Create
# a collection of set files which
# can be used later to traverse
# these sets for alignment.

import sys, os, shutil, random
from utils import scop, pfam, IO
from utils.logfile import XMLfile

if len(sys.argv) != 5: 
    print 'usage: %s SCOP_METADATA_PATH PDB_CACHE PDB_NUM_THRES SUPERFAMILY' % (sys.argv[0])
    exit()

scpath = sys.argv[1]
print '----> SCOP Metadata Path: %s' % (scpath)
pdcache = sys.argv[2]
print '----> PDB File Cache Path: %s' % (pdcache)
numthr = int(sys.argv[3])
print '----> Number of PDBs Threshold per Family: %d' % (numthr)
sfquery = sys.argv[4]
print '----> Superfamily Query: %s' % (sfquery)

# Get SCOP metadata.
scinst = scop.scopHierarchy(scpath)
print '====> Acquiring SCOP metadata <cache %s>' % (IO.getFolderName(scpath))
scinst.populateHierarchy()
print '====> Acquiring PDB files and set data'


# Acquire SCOP superfamily, if possible, from the given query.
print '====> Acquiring family data for query'
families = scinst.getFamilies(sfquery)
if len(families) == 0:
    print '----> No families were found given the query.'
    exit(1)
else:
    print '----> %d families were found given the query.' % (len(families))
    print '----> Recording list of families to %s.families.' % (sfquery)
o = open('%s.families' % (sfquery),'w')
for fam in families: o.write('%s\n' % (fam.sunid))
o.close()

# Iterate over families.
for fam in families:
    
    print '====> Acquiring PDBs for family <%s>' % (fam.shortname)
    print '----> <%s> has sunid %s' % (fam.shortname, fam.sunid)
    print '----> <%s> has description %s' % (fam.shortname, fam.desc)
    print '----> %d PDBs found...' % (len(fam.pdbs))
    if len(fam.pdbs) < numthr:
        print '----> This number of PDBs below threshold %d; ignoring family...' % (numthr)
        continue
    print '----> Downloading PDBs and saving to %s-%s.xml...' % (sfquery,fam.sunid)

    # Download necessary PDBs, first checking cache.
    xml = XMLfile('%s-%s.xml' % (sfquery,fam.sunid),'family')
    for pdb in fam.pdbs:
        o = pfam.grabPDBFile(pdb,pdcache)
        xml.add(xml.root,'pdb',('xml',o))
    