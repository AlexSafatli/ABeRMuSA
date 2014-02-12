# aligner.py
# ---------------------------------------------
# May 10, 2013; Alex Safatli
# ---------------------------------------------
# Deprecated threaded version of aligner.py.

import glob, os, cPickle
from threading import Thread
from utils import IO, PDBnet, homology, timer
from scipy.stats import norm

# Get scores from refscorer pickled file.

def getScores(args,ref,onlyAvg=False):
    reffile = open(os.path.join(IO.getFileName(ref),'ref.pickl'))
    scores, avgscr, highpv, rank = cPickle.load(reffile)
    reffile.close()
    if onlyAvg: return avgscr
    flatten = [scores[x] for x in scores]
    if len(flatten) != 0:
        minval,maxval = min(flatten),max(flatten)
        for key, val in scores.iteritems():
            if val == minval: mincol = key
            elif val == maxval: maxcol = key
    else: mincol,maxcol = -1,-1
    return (avgscr, mincol, maxcol, highpv, rank)

# Thread to manage execution of a command.

class executor(Thread):
    
    def __init__(self,args,ref,logf,exe,thrl):
        Thread.__init__(self)
        self.args = args
        self.ref = ref
        self.logf = logf
        self.exe = exe
        self.thrlist = thrl
    
    def run(self):
        
        # Add self to thread list.
        self.thrlist.append(self)

        # Execute the command.
        plugin = self.exe.plugin
        plugin.executeCmd(self.args,self.ref,self.exe,self.logf)

        # Remove self from THREAD_LIST.
        self.thrlist.remove(self)

# Thread to manage scoring of references.

class refscorer(Thread):
    
    def __init__(self,args,ref,logf,exe,thrl,thrd):
        Thread.__init__(self)
        self.args = args
        self.ref = ref
        self.log = logf
        self.exe = exe
        self.thrlist = thrl
        self.thrdone = thrd
        self.numHighPvals = 0

    def run(self):
        
        # Add self to thread list.
        self.thrlist.append(self)
        
        # Setup parameters.
        ref,args,exe = self.ref,self.args,self.exe
        reffldr = IO.getFileName(ref)
        files = '%s-%s.*.%s' % (self.exe.plugin.default_exe,
                                reffldr,self.exe.plugin.fasta_ext)
        filesl = glob.glob(os.path.join(reffldr,files))
        n = len(args)
        
        # While waiting for processing to complete, or even if complete, start
        # processing scores.
        failed = []
        scored = []
        scores = {}
        rmsddi = {}
        while len(scored) < (n-1):
            # Score a file. Adopted from MSA_pairwise_fester code.
            for fi in filesl:
                # Do not score if already scored.
                name = fi.strip(self.exe.plugin.fasta_ext).strip('.')
                if fi in scored: continue
                else:
                    try:
                        sc, rc = score(name + '.pdb',fi,self)
                    except IOError:
                        if not name in failed:
                            n -= 1
                            failed.append(name)
                        continue
                scored.append(fi)
                scores[name] = sc
                rmsddi[name] = rc
            filesl = glob.glob(os.path.join(reffldr,files))
        flatten = [scores[x] for x in scores if scores[x]]
        if len(flatten) == 0: avgscr = -1 # Insignificant score.
        else: avgscr = sum(flatten)/float(len(flatten))
        highpv = self.numHighPvals
        rank = avgscr*(highpv+1)
        
        # Write score vector to reference report.
        f = open('%s/ref.pickl' % reffldr, 'w')
        cPickle.dump((scores,avgscr,highpv,rank),f)
        f.close()
        
        # Write RMSD dictionary for analyses.
        f = open('%s/rmsds.pickl' % reffldr, 'w')
        cPickle.dump(rmsddi,f)
        f.close()
        
        # Report on success.
        self.log.write('%s had %d sequences (%s) that possessed sequence mismatches.' \
                       % (reffldr,len(failed),', '.join(failed)))
        self.log.write('Scored %s successfully; dumped to pickle file.' % (reffldr))
        self.log.incrementTimer()

        # Remove self from THREAD_LIST and add foldername to THREADS_DONE.
        self.thrlist.remove(self)
        self.thrdone.append(reffldr)

# Score an alignment.

def score(alignPDB,alignFASTA,refscr,scop=None):
    '''
    Gets alignment score irregardless of alignment method.
    '''
    
    # Get alignment length.
    alignlen = len(PDBnet.PDBstructure(alignPDB))
    
    # Get RRMSD and RMSD if length of alignment >= 100 residues.
    if not scop or alignlen >= 100:
        rrmsd, rmsd = homology.rrmsd(alignPDB,alignFASTA,True)
    
    # If non-significant P-value, return None.
    distr = norm(1,0.11)
    pval = distr.pdf(rrmsd)
    if pval > 0.05:
        refscr.numHighPvals += 1
        return None
    
    # Return the score value
    return rrmsd, rmsd

# Run an alignment.

def run(args,logf,ref=None,exe=None,quick=None,recursivecall=False):
    '''
    Runs an alignment; depends on a list of files to align, a logfile instance,
    an exewrapper instance, and has optional parameters.
    '''

    # Thread management: incl. list of currently-running threads,
    # a running queue of threads, a list of all threads done,
    # what the maximum number of threads in the thread list is,
    # and how many threads have launched.
    thrd_list, thrd_queue, thrd_done = [], [], []
    thrd_maxn, thrd_launched = 10, 0
    
    # Performance management.
    thrd_time = 0 # time since last thread done
    done_countd = 0 # number of threads done since last checked
    folders_out = []
    
    def handleThread():
        if len(thrd_list) < thrd_maxn and len(thrd_queue) > 0: 
            thrd_queue.pop(0).start()
            
    def handleFolder(ref):
        # Folder management.
        f = IO.getFileName(ref)
        if os.path.isdir(f):
            # Stop if already exists and contains items.
            return (not os.path.isfile('%s/ref.pickl' % (f)))
        else:
            # Make folder if it does not exist.
            os.mkdir(f)
            return True

    # Handle reference(s).
    if ref: refs = [ref] # A specific reference is provided.
    elif not quick: refs = args # No reference is provided; do n^2 search.
    else: refs = quick.get() # Use quick reference search method.

    # Create a sequence of threads to manage reference scoring and execute commands.
    for ref in refs:
        go = handleFolder(ref) # Folder management.
        if not go:
            folders_out.append(IO.getFileName(ref))
            continue
        # Execute command.
        r = refscorer(args,ref,logf,exe,thrd_list,thrd_done)
        e = executor(args,ref,logf,exe,thrd_list)
        thrd_queue.extend([r,e])
        thrd_launched += 1
    
    # Manage threads efficiently.
    while len(thrd_queue) > 0:
        # Wait until all threads complete.
        exe.cleanup() # cleanup any files.
        if (done_countd != 0 and done_countd < len(thrd_done)):
            if (timer.getTime() - thrd_time < 20) and (thrd_maxn < 100): 
                thrd_maxn += 1
                logf.write('Expanded thread queue by 1; thread count now %d.' \
                           % (thrd_maxn))
                thrd_time = timer.getTime()
        done_countd = len(thrd_done)
        handleThread()
    while len(thrd_done) < thrd_launched: pass
    
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
                        r = run(args,logf,ref=x,exe=exe,recursivecall=True)
                        refdict[x] = r[0][x]
                        thrd_done.extend(r[1])
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
    
    folders_out.extend(thrd_done)
    return (refdict, folders_out)

# Parse an alignment.

def parse(args,refdict,logf,mostsig=True):
    '''
    Parse an alignment run by the respective method.
    '''
    
    # See if anything to parse.
    if len(refdict) == 0:
        logf.write('Nothing to parse.')
        exit()
    
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
    logf.write('Best reference found to be %s (%s).' \
               % (IO.getFileName(bestref[0]),bestref[0]))
    logf.write('%s: rank measure %f.' \
               % (IO.getFileName(bestref[0]),bestref[1]))
    reffldr = IO.getFileName(bestref[0])

    return bestref, reffldr
