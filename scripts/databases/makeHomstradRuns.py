#!/usr/bin/env python

# makeHomstradRuns.py
# ------------------------------------
# Jul 5,  2013; Alex Safatli / Created
# Jan 13, 2013; Alex Safatli / Revamped
# ------------------------------------
# Walk through Homstrad Database
# and write Bash scripts in order
# to run ABeRMuSA and traditional
# alignments.

import sys, os
from utils import homstrad

NUM_QUICKSEARCH_ITERS = 20

if len(sys.argv) != 4: 
    print 'usage: %s MULTI HOMSTRAD_DB ALIGNER' % (sys.argv[0])
    exit()

ismulti = (sys.argv[1] == 'True')
print '----> Cluster Multi-Core Processing: %s' % (ismulti)
dbpath = sys.argv[2]
print '----> Database Path: %s' % (dbpath)
align = sys.argv[3]
print '----> Aligner Chosen: %s' % (align)
maxbin = 0

print '====> Reading Homstrad Database Files'
db = homstrad.homstradDatabase(dbpath)
db.traverse()
print 'Folders read: %d.' % (db.succeeded)
print '====> Writing Chain Run Commands to File'
fname = 'ABeRMuSA-%s-homstrad.sh' % (align)
shf = open(fname,'w')
if ismulti: m = ' -m 100 '
else: m = ' '

# ABeRMuSA commands.

for fl in db.folders:
    
    # Get all PDB paths.
    fl   = db.folders[fl]
    name = fl.name
    pdbs = [x for x in fl.getPDBs()]
    if (len(pdbs)) > maxbin: maxbin = len(pdbs)
    
    # Convert to long string; write to file as UNIX variable.
    pdbli = ' '.join(pdbs)
    pdbvn = 'PDB%s' % (name.upper())
    shf.write('%s="%s"\n' % (pdbvn,pdbli))
    
    # Ensure time elapsed is recorded to file.
    tp = '/usr/bin/time -o %s.elapsed' % (name)
   
    # Write commands.
    shf.write('%s ABeRMuSA -a %s -l%s-p %s $%s\n' % (tp,align,m,name,pdbvn))
    for i in xrange(1,NUM_QUICKSEARCH_ITERS+1):
        qprefix = name+'_q%d' % (i)
        shf.write('%s ABeRMuSA -a %s -q%d -l%s-p %s $%s\n' % (tp,align,i,m,qprefix,pdbvn))
        
shf.close()
print 'ABeRMuSA chain run commands saved to ./%s' % (fname)

# MSt (Multistructural) Alignment

# Matt
fn       = 'matt-msta-homstrad.sh'
shf      = open(fn,'w')
did, bid = {}, {}
bins     = range(0,maxbin+1)
for fl in db.folders:
    fl   = db.folders[fl]
    name = fl.name
    pdbs = [x for x in fl.getPDBs()]
    if (len(pdbs) in bins):
        bins.remove(len(pdbs))
        pdbli = ' '.join(pdbs)
        did[name] = pdbli
        bid[len(pdbs)] = name
        shf.write('echo \'---> running %s\'\n' % (name))
        shf.write('/usr/bin/time -o %s.elapsed Matt -o %s -t 2 %s > %s.stdout 2> %s.stdout\n' 
                  % (name,name,pdbli,name,name))
shf.close()
print 'Matt chain run commands saved to ./%s' % (fn)

# Write the number of structure --> group mappings.
fh = open('msta-homstrad.bins','w')
for num in bid: fh.write('%d\t%s\n' % (num,bid[num]))
fh.close()
print 'Number of structure to group mappings written to ./msta-homstrad.bins'