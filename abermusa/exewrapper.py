''' Encapsulates the execution of all pairwise alignments for a given reference structure. Handles generation of random jobnames if necessary for multiprocessing on a cluster. '''

# Date:   May 6 2013
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

import os, glob, subprocess
from labblouin.passToqsub import returnScript as qscript
from joblib import Parallel, delayed

# Function to check for existence of a binary on PATH.

def exeExists(cmd):
    return subprocess.call('type %s' % (cmd),shell=True,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE) == 0

# function ouside class for multiprocessor without a cluster
def individual_call(cmd):
    #print cmd # debugging
    sproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    err = sproc.communicate()
# Class definition(s).

class exewrapper:
    
    ''' Control centre for job execution and pairwise aligner behavior. '''
    
    def __init__(self,name,cmd,plugin,log,uniq=1,ismodel=False,ver='0.0.0'):
        
        ''' Setup control centre for job execution. '''
        
        self.uniq       = (uniq > 0)
        self.numcores   = uniq
        self.ismodel    = False
        self.name       = name
        self.cmd        = cmd
        self.plugin     = plugin
        self.logf       = log
        self.jobs       = []
        self.queue      = []
        self.ran        = {}
        self.fldr       = ''
        self.scpdbs     = None
        self.tag        = -1
        self.scoresToDo = None
        self.version    = ver
        
    def addScoreToDo(self,sco):
        
        ''' Add a score that should be computed for all alignments. '''
        
        if self.scoresToDo == None:
            self.scoresToDo = []
        self.scoresToDo.append(sco)
        
    def run(self):
        
        ''' Run all commands currently in the job queue. '''
            
        # Run all commands in the queue.
        cmds = [cmd[0] for cmd in self.queue]
        if not self.uniq:
            j = Parallel(n_jobs=-1)(delayed(individual_call)(cmd) for cmd in cmds)
        
        for cmd,fiout,fi,ref,o in self.queue:
            self.logf.incrementTimer()
            self.logf.writeTemporary('Aligned (%s, %s) to <%s>...' % (ref,o,fi))
            '''if self.uniq: cmds.append(cmd)
            else:
                sproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, shell=True)
                err = sproc.communicate()'''
            self.assertDone(ref,fiout,cmd)
            
        # Multiprocessing assumes a GRID Engine cluster.
        if self.uniq and len(cmds) > 0:
            if (not exeExists('qsub')):
                raise EnvironmentError('Multiprocessing assumes existence of qsub executable.')
            self.tag += 1
            jobname = '%s%d' % (self.name,self.tag)
            self.jobs.append(jobname)
            if len(cmds) == 1: cmd = qscript(cmds[-1],jobname)
            else: cmd = qscript(' ; '.join(cmds),jobname)
            sproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
            err = sproc.communicate()     
            
        # Clear the queue.
        self.queue = []
        
    def add(self,cmd,fiout,fi,ref,o):
        
        ''' Add a pairwise alignment to the job queue. '''
        
        self.queue.append((cmd,fiout,fi,ref,o))
        
    def assertDone(self,ref,fiout,cmd=''):
        
        ''' Forces the logic to treat a reference as done. '''
        
        if ref not in self.ran: self.ran[ref] = {}
        self.ran[ref][fiout] = cmd
        
    def cleanup(self):

        ''' Perform cleanup of all job/backup/GRID engine files. '''
        
        for job in self.jobs:
            toclean = glob.glob('*%s*' % (job))
            toclean.extend(glob.glob('core.*'))
            for fi in [x for x in toclean if os.path.isfile(x)]: os.remove(fi)