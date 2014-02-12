#!/usr/bin/env python

# mstaRMSD.py
# -------------------------
# Jul 22, 2013; Alex Safatli
# -------------------------
# Read PDBs from MSTA process
# and assign an RMSD to their
# alignments.

# Imports

import sys, os
from utils import IO, PDBnet, homology

# Constants

OUT_FILE = 'scores.csv'

# Script

if len(sys.argv) < 3:
    print 'usage: %s fasta_ext pdbfile [pdbfile2] ...' % (sys.argv[0])
    exit()

fasta_ext = sys.argv[1]
filist = sys.argv[2:]
rmsds = {}

for fi in filist:
    if not os.path.isfile(fi):
        print 'Not a file: <%s>. Skipping.' % (fi)
        continue
    fasta = IO.getFileName(fi) + '.%s' % (fasta_ext)
    if not os.path.isfile(fasta):
        print 'No FASTA file for: <%s>. Skipping.' % (fasta)
        rmsds[fi]  = -1
        continue
    try:
        p = PDBnet.PDBstructure(fi)
        rmsd  = p.rmsd(fasta)
    except Exception,e:
        print 'Could not score <%s>; reason: %s' % (fi,e)
        rmsds[fi]  = -1   
        continue
    print 'RMSD for <%s> is %f.' % (fi, rmsd)
    rmsds[fi]  = rmsd
    
o = open(OUT_FILE,'w')
o.write('file,rmsd\n')
for key in rmsds: o.write('%s,%f\n' % (key,rmsds[key]))
o.close()
print 'Results written to %s.' % (OUT_FILE)
    
    