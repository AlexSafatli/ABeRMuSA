# exewrapper.py
# ---------------------------------------------
# May 6, 2013; Alex Safatli
# ---------------------------------------------
# Encapsulates the execution of all pairwise alignments
# for a given reference structure. Handles generation
# of random jobnames if necessary for Fester multiprocessing.

import os, glob, subprocess
from utils.passToqsub import returnScript as qscript

class exewrapper:
    def __init__(self,name,cmd,plugin,log,uniq=1):
        self.uniq = (uniq > 0)
        self.numcores = uniq
        self.name = name
        self.cmd = cmd
        self.plugin = plugin
        self.logf = log
        self.jobs = []
        self.queue = []
        self.ran = {}
        self.fldr = ''
        self.scpdbs = None
        self.tag = -1
    def run(self):
        # Run all commands in the queue.
        cmds = []
        for cmd,fiout,fi,ref,o in self.queue:
            self.logf.incrementTimer()
            self.logf.write('Aligning (%s, %s) to <%s>...' % (ref,o,fi))
            if self.uniq: cmds.append(cmd)
            else:
                sproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, shell=True)
                err = sproc.communicate()
            self.assertDone(ref,fiout,cmd)
        if self.uniq and len(cmds) > 0:
            # assumes fester qscript support
            self.tag += 1
            jobname = '%s%d' % (self.name,self.tag)
            self.jobs.append(jobname)
            if len(cmds) == 1: cmd = qscript(cmds[-1],jobname)
            else: cmd = qscript(' ; '.join(cmds),jobname)
            sproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
            err = sproc.communicate()     
        self.queue = []
    def add(self,cmd,fiout,fi,ref,o):
        self.queue.append((cmd,fiout,fi,ref,o))
    def assertDone(self,ref,fiout,cmd=''):
        if ref not in self.ran: self.ran[ref] = {}
        self.ran[ref][fiout] = cmd
    def cleanup(self):
        for job in self.jobs:
            toclean = glob.glob('*%s*' % (job))
            toclean.extend(glob.glob('core.*'))
            for fi in [x for x in toclean 
                       if os.path.isfile(x)]: os.remove(fi)