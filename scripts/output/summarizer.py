#!/usr/bin/env python

# summarizer.py
# -------------------------
# Jul 22, 2013; Alex Safatli
# -------------------------
# Analyze all XML file output
# from Damastes in order to
# summarize results.

# Importing ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

import sys, os, pickle, numpy, cPickle
from utils import PDBnet, logfile, IO
from utils.logfile import parseTime

# Constants ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

OUT_FILE    = 'summary'
FAILED_FILE = 'failed'
DELETE_FILE = 'delete.sh'

# Damastes Run :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class damastesRun:
    def __init__(self,aXML):
        self.xml = aXML
        self.files = []
        self.reference = None
        self.reffldr = ''
        self.__read__()
    def __iter__(self):
        for fi in self.files: yield fi
    def __read__(self):
        xml = logfile.XMLfile(self.xml)
        xml.read()
        if xml.root.tag != 'alignment':
            raise IOError('Invalid ABeRMuSA XML file.')
        for child in xml.root:
            if child.tag == 'file': self.files.append(child.text)
            elif child.tag == 'reference': self.reference = child.text
        if self.reference: self.reffldr = IO.getFileName(self.reference)
    def getReferences(self,damastesLog):
        # Open a Damastes logfile and
        # extract references in order
        # they appear.
        refs = []
        if not os.path.isfile(damastesLog): return []
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
    def getElapsed(self,damastesLog):
        # Open a Damastes logfile, get
        # elapsed time.
        if not os.path.isfile(damastesLog): return "-"
        f = open(damastesLog)
        li = f.readlines()
        f.close()
        time = "-"
        for l in li:
            if '] Elapsed:' in l:
                find = l.find('Elapsed:')
                ttime = l[find+9:].strip()
                if time == "-" or parseTime(ttime) > parseTime(time):
                    time = ttime
        return time

# Script Body ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print 'usage: %s [-includeq] [damastesXMLfile1] [damastesXMLfile2] ...' % (sys.argv[0])
        exit()
    
    # Assumes naming convention where if _q# is in
    # filename of XML file, it is a quick-search mode
    # run.
    
    incflag = False
    runs   = {}
    scores = {}
    
    # Gather all data.
    
    for xml in sys.argv[1:]:
        # For every XML file given from Damastes...
        if xml.lower() == '-includeq':
            incflag = True
            print 'Recognized flag for including quick search XMLs.'
            continue
        run = damastesRun(xml)
        if not '_q' in xml: 
            xmlind = xml.find('.xml')
            key = xml[:xmlind]
            if not key in runs: runs[key] = {}
            runs[key]['exh'] = run
        elif incflag:
            xmlind = xml.find('.xml')
            key = xml[:xml.find('_q')]
            if not key in runs: runs[key] = {}
            runs[key][xml[xml.find('_q')+2:xmlind]] = run            
        else: continue
    
    # Report results.
    
    failed = []
    fh = open(OUT_FILE,'w')
    fh.write('Group\tNum. Structs\tElapsed\tOutput\tRMSD\n')
    numstructs = 0
    for r in runs:
        if 'exh' in runs[r]: arun = runs[r]['exh']
        else: # Get highest possible iteration count.
            keys = [int(x) for x in runs[r].keys()]
            highestkey = str(max(keys))
            arun = runs[r][highestkey]
        numstructs += len(arun.files)
        hasoutput = 'Y'
        if (arun.reference == None):
            failed.append(r)
            hasoutput = 'N'
            avgrmsd = "-"
        else:
            reffile = os.path.join(arun.reffldr,'ref.pickl')
            o = open(reffile)
            scores, _, _, _ = cPickle.load(o)
            o.close()
            armsds = []
            for key in [x for x in scores if scores[x]]:
                scfi = '%s/%s.%s' % (arun.reffldr,key,'pickl')
                o = open(scfi)
                _, rc, _ = cPickle.load(o)
                o.close()
                armsds.append(rc)
            avgrmsd = "%f" % numpy.average(armsds)
        logfn = IO.getFileName(r)
        name = r
        if not 'exh' in runs[r]:
            logfn += '_q%s' % (highestkey)
            name +=  '_q%s' % (highestkey)
        logfn += '.log'
        fh.write('%s\t%d\t%s\t%s\t%s\n' % (name,len(arun.files),
                                           arun.getElapsed(logfn),
                                           hasoutput, avgrmsd))
        
    if len(failed) > 0:
        df = open(DELETE_FILE,'w')
        fh = open(FAILED_FILE,'w')
        for fail in failed:
            arun = runs[fail][runs[fail].keys()[0]]
            todelete = [IO.getFileName(x) for x in arun.files]
            for d in todelete: df.write('rm %s/*.pickl 2> /dev/null\n' % (d))
            fh.write('%s\n' % (fail))
        fh.close()
        df.close()
        
    print '\nDataset size: %d' % (len(runs))
    print 'Number of structures: %d' % (numstructs)
    print 'Number of failed runs: %d' % (len(failed))
    print 'DONE.'    
    fh.close()
    




