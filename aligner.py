''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++

aligner.py
May 10th, 2013; Alex Safatli

Carries out alignment and scoring. Fester/multiprocessor execution is abstracted through the exewrapper class (see appropriate file). Running an alignment transitions between two states: execution of commands and scoring of output. '''

import glob, os, cPickle, math, scoring
from utils import IO, PDBnet, homology, timer
from utils.logfile import XMLfile
from scipy.stats.kde import gaussian_kde

# Get scores from pickled file.

def getScores(args,ref,onlyAvg=False):
    
    ''' Open a reference pickled file and retrieve information. '''   
    
    reffile = open(os.path.join(IO.getFileName(ref),'ref.pickl'))
    scores, avgscr, highpv, rank, sctype = cPickle.load(reffile)
    reffile.close()
    if onlyAvg: return avgscr
    mincol,maxcol = -1,-1
    flatten = [scores[x] for x in scores]
    if sctype != 'RMSD' and sctype != 'RRMSD':
        flatten = [(1-x) for x in flatten] # Flip directions of floats.
    if len(flatten) != 0:
        minval,maxval = min(flatten),max(flatten)
        for key, val in scores.iteritems():
            if val == minval:   mincol = key
            elif val == maxval: maxcol = key
    return (avgscr, mincol, maxcol, highpv, rank)

def handleFolder(f):
    
    ''' Determine if a reference folder has already been scored. '''
    
    if os.path.isdir(f):
        # Stop if already exists and contains items.
        return (not os.path.isfile('%s/ref.pickl' % (f)))
    else:
        # Make folder if it does not exist.
        os.mkdir(f)
        return True

def handleReference(args,ref,exe,logf):
    
    ''' Score a reference folder. '''
    
    # Preparation steps.
    rf = IO.getFileName(ref)
    if not handleFolder(rf): return False
    elif rf not in exe.ran: return False
    fl = exe.ran[rf].keys() # List of files.
    
    # Data structure initialization.
    failed,fdetails,scored,scores,pvals = [],[],[],{},{}
    sctype = None
    
    # For all files...
    for f in fl:
        
        # Get metadata.
        name = IO.getFileName(f)
        pr   = '%s/%s' % (rf,name) # Prefix.
        scfi = '%s.pickl' % (pr)   # Score file.
        pdbf = '%s.pdb'   % (pr)   # PDB file.
        outf = '%s.out'   % (pr)   # Output file from plugin aligner.
        
        # See if done yet.
        if os.path.isfile(scfi): scored.append((name,scfi)) # Scored already.
        elif not os.path.isfile(outf) and not os.path.isfile(f): continue # Not done yet.
        elif not os.path.isfile(f):
            # Output was generated by program; probably failed.
            failed.append(name)
            fdetails.append('Failed scoring %s. Reason: missing <%s>.' % (name,f))
        elif not os.path.isfile(scfi) and not name in failed:
            try:
                # Post-process (plugin-dependent) and score.
                if exe.plugin.postprocess: exe.plugin.postProcess(pdbf,ref,logf)
                sc = scoring.score(pdbf,f,exe,logf)
                if not sctype: sctype = sc[0].getName()
                scores[name] = sc[0].getScore()
                pvals[name]  = sc[0].getPValue()
                scoref = scoring.scoreFile(scfi)
                scoref.addScores(sc)
                scoref.writeToFile()
            except Exception,e:
                fdetails.append('Failed scoring %s. Reason: %s.' % (name,str(e)))
                failed.append(name)
                continue
            logf.writeTemporary('Scored %s <%s: %f, p-value: %s>.' % (
                name,sctype.lower(),scores[name],str(pvals[name])))
            try: os.remove(outf) # Removes its stdout file if succeeded.
            except: pass
            scored.append((name,scfi))
            logf.incrementTimer()
    
    # Check to see if done all needed files.
    if (len(scored) + len(failed)) >= (len(args)-1):
        
        # Report all failures (this is to ensure it is only done once).
        logf.write('%s had %d sequences (%s) that failed to align.' 
                   % (rf,len(failed),', '.join(failed)))        
        for fd in fdetails: logf.write(fd)
        
        # Score all pairwise alignments together for the reference.
        highpv = 0 # Number of high p-values.
        for name,scfi in scored:
            # If not already cached...
            if not name in scores:
                scoref = scoring.scoreFile(scfi)
                sc = scoref.getScores()
                if not sctype: sctype = sc[0].getName()
                scores[name] = sc[0].getScore()
                pvals[name]  = sc[0].getPValue()
            # See if low P-value.
            if pvals[name] and pvals[name] < 0.05:
                highpv += 1
                scores[name] = None
                
        # Flatten all RRMSD values into a linear array.
        flatten = [scores[x] for x in scores if scores[x]]
        if sctype != 'RMSD' and sctype != 'RRMSD':
            flatten = [(1-x) for x in flatten] # Flip directions of floats.
        if len(flatten) == 0: avgscr = -1
        else: avgscr = sum(flatten)/float(len(flatten))
        rank = avgscr*(highpv+1) # Calculate rank.
        
        # Write score vector reference report.
        picklf = '%s/ref.pickl' % (rf)
        o = open(picklf,'w')
        cPickle.dump((scores,avgscr,highpv,rank,sctype),o)
        o.close()
        logf.write('Coverage of reference complete (%d succeeded, %d failed).' 
                   % (len(scored),len(failed)))
        
        # Write manifest for this reference.
        xmlf = '%s/manifest.xml' % (rf)
        xml = XMLfile(xmlf,'reference')
        xml.add(xml.root,'score',('type',sctype))
        for scoredn,_ in scored:
            xml.add(xml.root,'succeeded',('xml',scoredn),('score',scores[scoredn]))
        for failedn in failed:
            xml.add(xml.root,'failed',   ('xml',failedn))
        logf.write('Manifest XML file written for %s; located at <%s>.' % (rf,xmlf))
        
        # Report success.
        logf.write('Scored %s successfully; dumped to pickle file <%s>.' % (rf,picklf))
        return True
    
    return False

# Run an alignment.

def run(args,logf,ref=None,exe=None,quick=None,recursivecall=False):
    
    ''' Runs an alignment; depends on a list of files to align, a logfile instance,
    an exewrapper instance, and has optional parameters. '''

    # Determine what is done.
    folders_out, i = [], 0

    # Handle reference(s).  
    if ref: refs = [ref]          # A specific reference is provided.
    elif not quick: refs = args   # No reference is provided; do n^2 search.
    else: refs = quick.get()      # Use quick reference search method.

    def checkAll():
        for r in [x for x in refs if not IO.getFileName(x) in folders_out]:
            f = handleReference(args,r,exe,logf)
            if f: folders_out.append(IO.getFileName(r))
        exe.cleanup() # cleanup any files

    # Report reference list.
    for ref in refs: logf.write('Reference: %s' % (ref))

    # Execute and score for every reference.
    for ref in refs:
        i += 1
        reffldr = IO.getFileName(ref)
        go = handleFolder(reffldr) # Is already done?
        if not go:
            folders_out.append(reffldr)
            logf.writeTemporary('Reference (%s) alignment already done.' % (reffldr))
            logf.updateTimer(logf.numat+2*(len(args)-1))
            continue
        # Execute command; abstracts multi-processing if necessary.
        exe.plugin.executeCmd(args,ref,exe,logf)
        if ((exe.numcores == 0) or ((i % exe.numcores) == 0)): checkAll()
    
    # If done executing (e.g., all to Fester), wait around and score as files
    # are output from running alignments.
    while (len(folders_out) != len(refs)): checkAll()
    
    # Be prepared to store all average scores for references.
    refdict = {}
    
    # If quick reference, have to manage further executions.
    if quick:
        for ref in refs:
            # Get score vector from file.
            avgscr, minref, maxref, _, _ = getScores(args,ref)
            if avgscr >= 0:
                refdict[ref] = avgscr
                for x in [minref,maxref]:
                    if x != -1: # If x == -1, nothing was found.
                        fargs = [IO.getFileName(i) for i in args]
                        toref = x[x.rfind('.')+1:]
                        if toref in fargs: recref = fargs.index(toref)
                        else:
                            logf.write('WARNING; Quick search reference <%s> missing.' % (toref))
                            continue
                        if args[recref] in refdict: continue
                        r = run(args,logf,ref=args[recref],exe=exe,recursivecall=True)
                        if len(r[0]) > 0:
                            refdict[args[recref]] = r[0][args[recref]]
                            folders_out.extend(r[1])
            else: logf.write('%s found to be insignificant.' % (IO.getFileName(ref)))
    else:
        for ref in refs:
            avgscr = getScores(args,ref,True)
            if avgscr >= 0: refdict[ref] = avgscr
            else: logf.write('%s found to be insignificant.' % (IO.getFileName(ref)))             
    
    if not recursivecall:
        print '' # Newline to console.
        for ref in refdict:
            # Report on findings to log file.
            logf.write('%s average score found to be %f.' \
                       % (IO.getFileName(ref),refdict[ref]))
            
    return (refdict, folders_out)

# Parse an alignment.

def parse(args,refdict,logf,mostsig=True):
   
    ''' Parse an alignment run by the respective method. '''
    
    # See if anything to parse.
    if len(refdict) == 0:
        logf.write('ERROR; No significant references found. Nothing to parse!')
        return None, None
    
    # Determine the best reference.
    bestref = None
    for ref in refdict:
        avg, _, _, highpv, rank = getScores(args,ref)
        if not bestref: bestref = [ref, rank, highpv]
        else:
            if mostsig and (highpv < bestref[2] or \
                            (highpv == bestref[2] \
                             and rank < bestref[1])):
                    bestref = [ref, rank, highpv]
            elif not mostsig and rank < bestref[1]:
                bestref = [ref, rank, highpv]
    reffldr = IO.getFileName(bestref[0])
    logf.write('Best reference found to be %s (%s).' % (reffldr,bestref[0]))
    logf.write('%s: rank measure %f.\n' % (reffldr,bestref[1]))
    
    return bestref, reffldr
