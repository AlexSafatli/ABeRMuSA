#!/usr/bin/env python

# quickModeTester.py
# -------------------------
# Jul 5, 2013; Alex Safatli
# -------------------------
# Analyze all XML file output
# from ABeRMuSA in order to
# analyze the performance of the
# quick search heuristic/mode.

# Importing ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

import sys, os, pickle, numpy, cPickle
try:
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt
    plotEnabled = True
except Exception,e:
    print 'WARNING; library "matplotlib" not found. Plotting disabled.'
    plotEnabled = False
from utils import PDBnet, logfile, IO

# Constants ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

PICKLE_FILE = 'quickModeTester.pickl'

# Damastes Run :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class damastesRun:
    def __init__(self,damastesXML):
        self.xml = damastesXML
        self.files = []
        self.reference = None
        self.reffldr = ''
        self.__read__()
    def __read__(self):
        xml = logfile.XMLfile(self.xml)
        xml.read()
        if xml.root.tag != 'alignment': raise IOError('Invalid ABeRMuSA XML file.')
        for child in xml.root:
            if child.tag == 'file': self.files.append(child.text)
            elif child.tag == 'reference': self.reference = child.text
        if self.reference: self.reffldr = IO.getFileName(self.reference)
    def getReferences(self,damastesLog):
        # Open a Damastes logfile and
        # extract references in order
        # they appear.
        refs = []
        f = open(damastesLog)
        li = f.readlines()
        f.close()
        for l in li:
            if '] Reference:' in l:
                find = l.find('Reference:')
                ref  = l[find+10:].strip()
                refs.append(ref)
            elif 'WARNING; Reference <' in l:
                find = l.find('<')
                ref  = l[find+1:l.find('>')].strip()
                refs.append(ref)
        return refs

# Reference Tree :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class refNode:
    def __init__(self,ref):
        self.ref = ref
        self.note = ''
        self.next = []
    def addNext(self,n):
        self.next.append(n)
    def __iter__(self):
        for n in self.next: yield n

class refTree:
    def __init__(self,refs,iters):
        self.refs = refs
        self.iters = iters
        self.root = None
        self.__build__()
    def __build__(self):
        self.root = refNode(None)
        for r in range(0,self.iters):
            ref = refNode(self.refs[r])
            ref.note = 'initial'
            self.root.addNext(ref)
            for n in range(self.iters+r,self.iters+r+2):
                if n in range(0,len(self.refs)):
                    s = refNode(self.refs[n])
                    if (self.iters+r+3-n == 0): s.note = 'min'
                    else: s.note = 'max'
                    ref.addNext(s)

# Script Body ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if len(sys.argv) < 2:
    print 'usage: %s [damastesXMLfile1] [damastesXMLfile2] ...' % (sys.argv[0])
    exit()

# Assumes naming convention where if _q# is in
# filename of XML file, it is a quick-search mode
# run.

runs = {}
found_in = []
found_as = []
scores   = {}
never_fd = 0

# Gather all data.

for xml in sys.argv[1:]:
    # For every XML file given from Damastes...
    run = damastesRun(xml)
    if not '_q' in xml: 
        xmlind = xml.find('.xml')
        key = xml[:xmlind]
        if not key in runs: runs[key] = {}
        runs[key]['exh'] = run
    else:
        xmlind = xml.find('.xml')
        key = xml[:xml.find('_q')]
        if not key in runs: runs[key] = {}
        runs[key][xml[xml.find('_q')+2:xmlind]] = run

# Pickle this dictionary.
    
print '\nSaving found results to disk <%s>.' % (PICKLE_FILE)
fh = open(PICKLE_FILE,'w')
pickle.dump(runs,fh)
fh.close()

# Determine quick search performance.

print '\nDetermining quick search performance...'
for key in runs:
    # See how finding true reference worked.
    got = False
    run = runs[key]
    if 'exh' in run: bestref = run['exh'].reference
    else:
        print 'WARNING; No exhaustive search found for %s.' % (key)
        continue
    if not bestref: continue
    refpickl = open('%s/ref.pickl' % (run['exh'].reffldr))
    _, avg, _, _ = cPickle.load(refpickl)
    refpickl.close()
    bestrefavg = avg
    qs = sorted([int(x) for x in run if x != 'exh'])
    refcache = {}
    for q in qs:
        if not got and run[str(q)].reference == bestref:
            try: found_in.append(q)
            except: print 'WARNING; %s could not be parsed to integer.' % (str(q))
            got = True
            # See how reference picking behave.
            refs = run[str(q)].getReferences('%s_q%s.log' % (key,str(q)))
            rtre = None
            try: rtre = refTree(refs,q)
            except Exception,e:
                print 'WARNING; Could not find index %s in %s (%s).' \
                      % (str(q),refs,'%s_q%s'%(key,str(q)))
            if rtre:
                for r in rtre.root:
                    if r.ref == bestref: found_as.append(r.note)
                    else:
                        for x in r:
                            if x.ref == bestref: found_as.append(x.note)
        if run[str(q)].reference: 
            reftoget = IO.getFileName(run[str(q)].reference)
            if reftoget in refcache: avg = refcache[reftoget]
            else:
                refpickl = open('%s/ref.pickl' % (reftoget))
                _, avg, _, _ = cPickle.load(refpickl)
                refpickl.close()
                refcache[reftoget] = avg
            if not str(q) in scores: scores[str(q)] = []
            scores[str(q)].append(avg-bestrefavg)
    if not got: never_fd += 1

# Plot results.

if plotEnabled:
    plt.hist(found_in)
    print '\nPlotted and saved plots to disk <%s>.' % ('found_in.pdf')
    pp = PdfPages('found_in.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()

# Report results.

fh = open('found_in.txt','w')
fa = open('found_as.txt','w')
print '' # Newline.
perc = 0 # Reset percentage.
for val in set(found_in):
    count = found_in.count(val)
    perc  += (float(count)/len(found_in))*100
    avg   = numpy.average(scores[str(val)])
    print 'True reference found in iteration [n=%d] %d times (%d%%); average delta RRMSD: %f' \
          % (val,count,perc,avg)
    fh.write('%d\t%d\t%f\n' % (val,perc,avg))
print 'Length of testset: %d' % (len(found_in))
print '' # Newline.
for note in set(found_as):
    count = found_as.count(note)
    perc  = (float(count)/len(found_as))*100
    print 'True reference found in %s along quick search track %d times (%d%%).' \
          % (note,count,perc)
    fa.write('%s\t%d\t%f\n' % (note,count,perc))
print 'Length of testset: %d' % (len(found_as))

numstructs = 0
for r in runs:
    arun = runs[r][runs[r].keys()[0]]
    numstructs += len(arun.files)
print '\nDataset size: %d' % (len(runs))
print 'Number of structures: %d' % (numstructs)
print '%d groups never had a true reference found by the heuristic.' % (never_fd)
print '\nDONE.'    
fh.close()
fa.close()
    




