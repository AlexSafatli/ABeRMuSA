#!/usr/bin/env python

# makeSabmarkRuns.py
# -------------------------
# Jun 14, 2012; Alex Safatli
# -------------------------
# Walk through Sabmark Database
# and create Bash files for alignments.

import sys, os
from utils import sabmark

NUM_QUICKSEARCH_ITERS = 20

if len(sys.argv) != 5: 
    print 'usage: %s MULTI SABMARK_DB SABMARK_CAT ALIGNER' % (sys.argv[0])
    exit()

ismulti = (sys.argv[1] == 'True')
print '----> Cluster Multi-Core Processing: %s' % (ismulti)
dbpath  = sys.argv[2]
print '----> Database Path: %s' % (dbpath)
dbcat   = sys.argv[3]
print '----> Database Category: %s' % (dbcat)
align   = sys.argv[4]
print '----> Aligner Chosen: %s' % (align)
maxbin  = 0

print '====> Reading SABMARK Database Files'
db = sabmark.sabmarkDatabase(dbpath)
db.traverse()
print 'Folders read: %d' % (len(db.folders))
print '====> Writing Damastes Chain Run Commands to File'
print 'PDB records and validity saved to ./%s.db' % (dbcat)
print 'Damastes chain run commands to be saved to ./ABeRMuSA-%s.sh' % (dbcat)
print 'Matt MStA chain run commands to be saved to ./matt-msta-%s.sh' % (dbcat)
pdf = open('%s.db' % (dbcat),'w')
shf = open('ABeRMuSA-%s.sh' % (dbcat),'w')
if ismulti: m = ' -m 100 '
else: m = ' '

# Damastes :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for gr in db.groups((os.path.join(dbpath,dbcat),dbcat)):
    name = gr.name
    pdbs = gr.pdbs
    if (len(pdbs)) > maxbin: maxbin = len(pdbs)
    pdbli = ' '.join([x.path for x in pdbs])
    for pdb in pdbs: pdf.write('%s\t%s\n' % (pdb.path,pdb.true))

    # Convert to long string; write to file as UNIX variable.
    pdbli = ' '.join(pdbs)
    pdbvn = 'PDB%s' % (name.upper())
    shf.write('%s="%s"\n' % (pdbvn,pdbli))
    
    # Ensure time elapsed is recorded to file.
    if not ismulti: tp = '/usr/bin/time -o %s.elapsed' % (name)
    else:           tp = ''
   
    # Write commands.
    shf.write('%s ABeRMuSA -a %s -l%s-p %s $%s\n' % (tp,align,m,name,pdbvn))
    for i in xrange(1,NUM_QUICKSEARCH_ITERS+1):
        qprefix = name+'_q%d' % (i)
        shf.write('%s ABeRMuSA -a %s -q%d -l%s-p %s $%s\n' % (tp,align,i,m,qprefix,pdbvn))

shf.close()
pdf.close()

# MSt (Multistructural) Alignment ::::::::::::::::::::::::::::::::::

did = {}
bid = {}
bins = range(0,maxbin+1)
# Matt
shf = open('matt-msta-%s.sh' % (dbcat),'w')
for gr in db.groups((os.path.join(dbpath,dbcat),dbcat)):
    name = gr.name
    pdbs = gr.pdbs
    pdbli = ' '.join([x.path for x in pdbs])
    if len(pdbs) in bins and not os.path.isfile('%s.txt' % (name)):
        bins.remove(len(pdbs))
        did[name] = pdbli
        bid[len(pdbs)] = name
        shf.write('echo \'----> running %s\'\n' % (name))
        if not ismulti: timeprefix = '/usr/bin/time -o %s.elapsed ' % (name)
        else: timeprefix = ''
        shf.write('%sMatt %s -t 2 %s > %s.stdout 2> %s.stderr\n' 
                  % (timeprefix,name,pdbli,name,name))
shf.close()

# Write the number of structure --> group mappings.
fh = open('msta-%s.bins' % (dbcat),'w')
for num in bid: fh.write('%d\t%s\n' % (num,bid[num]))
fh.close()
print 'Number of structure to group mappings written to ./msta-%s.bins' % (dbcat)