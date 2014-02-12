#!/usr/bin/env python

# makeSCOPRun.py
# -------------------------
# Jul 5, 2013; Alex Safatli
# -------------------------
# Read XML SCOP set
# and wrap as shell script files
# to perform alignments.

import sys, os
from utils import IO
from utils.logfile import XMLfile

if len(sys.argv) < 4: 
    print 'usage: %s MULTI ALIGNER [file.xml] [file2.xml]' % (sys.argv[0])
    exit()

ismulti = (sys.argv[1] == 'True')
print '----> Cluster Multi-Core Processing: %s' % (ismulti)
align = sys.argv[2]
print '----> Aligner Chosen: %s' % (align)
xmlfiles = []
print '====> Reading Input'
for item in sys.argv[3:]:
    if not os.path.isfile(item):
        print 'Not an item. Skipping.'
        continue
    xml = XMLfile(item)
    xml.read()
    if xml.root.tag not in ['superfamily','family']: continue
    xmlfiles.append(xml)
maxbin_d = 0
maxbin_s = 0

print '====> Writing Chain Run Commands to File'
shf = open('ABeRMuSA-%s-scop.sh' % (align),'w')
if ismulti: m = ' -m 100 '
else: m = ' '

# Damastes :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for xml in xmlfiles:
    name = IO.getFileName(xml.path)
    name = ''.join(name.split(' '))
    pdbs_u = [] # Unclassified
    pdbs_d = [] # Dissimilar
    pdbs_s = [] # Similar
    for child in xml.root:
        if child.tag == 'dissimilar': pdbs_d.append(child.text)
        elif child.tag == 'similar': pdbs_s.append(child.text)
        else: pdbs_u.append(child.text)
    if (len(pdbs_d)) > maxbin_d: maxbin_d = len(pdbs_d)
    if (len(pdbs_s)) > maxbin_s: maxbin_s = len(pdbs_s)
    pdbli_d = ' '.join(pdbs_d)
    pdbli_s = ' '.join(pdbs_s)
    pdbli_u = ' '.join(pdbs_u)
    shf.write('echo \'----> running %s\'\n' % (name))
    if len(pdbs_d) > 0:
        shf.write('ABeRMuSA -x -a %s -l%s-p %s %s\n' % 
                  (align,m,'%s-%s_d' % (name,align),pdbli_d))
    if len(pdbs_s) > 0:
        shf.write('ABeRMuSA -x -a %s -l%s-p %s %s\n' % 
                  (align,m,'%s-%s_s' % (name,align),pdbli_s))
    if len(pdbs_u) > 0:
        shf.write('ABeRMuSA -x -a %s -l%s-p %s %s\n' % 
                  (align,m,'%s-%s' % (name,align),pdbli_u))
    for i in xrange(1,7):
        if len(pdbs_d) > 0:
            shf.write('ABeRMuSA -x -a %s -q %d -l%s-p %s %s\n' 
                      % (align,i,m,name+'-%s_d_q%d'% (align,i),pdbli_d))
        if len(pdbs_s) > 0:
            shf.write('ABeRMuSA -x -a %s -q %d -l%s-p %s %s\n' 
                      % (align,i,m,name+'-%s_s_q%d'% (align,i),pdbli_s))
        if len(pdbs_u) > 0:
            shf.write('ABeRMuSA -x -a %s -q %d -l%s-p %s %s\n' 
                      % (align,i,m,name+'-%s_q%d'% (align,i),pdbli_u))
shf.close()
print 'ABeRMuSA chain run commands saved to ./ABeRMuSA-%s-scop.sh' % (align)

# MSt (Multistructural) Alignment :::::::::::::::::::::::::::::::::::::::::::::::

# Matt
fn = 'matt-msta-scop.sh'
shf = open(fn,'w')
did_d = {}
bid_d = {}
bins_d = range(0,maxbin_d+1)
did_s = {}
bid_s = {}
bins_s = range(0,maxbin_s+1)
for xml in xmlfiles:
    name = IO.getFileName(xml.path)
    pdbs = []
    pdbs_d = []
    pdbs_s = []
    for child in xml.root:
        if child.tag == 'dissimilar': pdbs_d.append(child.text)
        elif child.tag == 'similar': pdbs_s.append(child.text)
    if (len(pdbs_d) in bins_d):
        bins_d.remove(len(pdbs_d))
        pdbli = ' '.join(pdbs_d)
        did_d[name] = pdbli
        bid_d[len(pdbs_d)] = name
        shf.write('echo \'---> running %s\'\n' % (name))
        shf.write('/usr/bin/time -o %s_d.elapsed Matt -o %s_d -t 2 %s > %s_d.stdout 2> %s_d.stdout\n' 
                  % (name,name,pdbli,name,name))
    if (len(pdbs_s) in bins_s):
        bins_s.remove(len(pdbs_s))
        pdbli = ' '.join(pdbs_s)
        did_s[name] = pdbli
        bid_s[len(pdbs_s)] = name
        shf.write('echo \'---> running %s\'\n' % (name))
        shf.write('/usr/bin/time -o %s_s.elapsed Matt -o %s_s -t 2 %s > %s_s.stdout 2> %s_s.stdout\n' 
                  % (name,name,pdbli,name,name))
shf.close()
print 'Matt chain run commands saved to ./%s' % (fn)

# Write the number of structure --> group mappings.
fh = open('msta-scop_d.bins','w')
for num in bid_d: fh.write('%d\t%s\n' % (num,bid_d[num]))
fh.close()
fh = open('msta-scop_s.bins','w')
for num in bid_s: fh.write('%d\t%s\n' % (num,bid_s[num]))
fh.close()
print 'Number of structure to group mappings written to ./msta-scop_X.bins'