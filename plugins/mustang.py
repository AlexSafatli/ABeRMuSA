# mustang.py
# ---------------------------------------------
# May 10, 2013; Alex Safatli
# ---------------------------------------------
# Plugin for ABeRMuSA adding support for the 
# Mustang pairwise aligner executable.

import os
from utils import IO

postprocess = False
default_exe = 'mustang'
fasta_ext   = 'afasta'

'''   
def scoreAlignment(refscr,pre):
    
    f = open(pre + '.rms_rot')
    line = f.readline()
    while not line.startswith('  1'):
        line = f.readline()
    rmsd = float(line.split()[-1])
    f.close()
    # Get length of alignment
    lenAlign = 0
    f = open(pre + fasta_ext)
    chains = ['', '']
    line = f.readline()
    line = f.readline()
    while not line.startswith('>'):
        chains[0] += line.strip()
        line = f.readline()
    for line in f:
        chains[1] += line.strip()
    if len(chains[0]) != len(chains[1]):
        print 'Length of the two chains in ' + prefix + ' are not the same.'
        exit()
    for q in range(len(chains[0])):
        if chains[0][q] != '-' and chains[1][q] != '-':
            lenAlign += 1
    if lenAlign == 0:
        return None
    # Return the score value
    return rmsd/float(lenAlign)
'''

def postProcess(fi,ref,log):
    
    # Nothing to do.
    
    pass

def executeCmd(args,ref,exe,log):
    
    # Decide on a reference folder.
    
    reffldr = IO.getFileName(ref)
    
    # Loop and add commands to exewrapper.
    
    for fi in [x for x in args if x != ref]:
        # For every file except the reference.
        outpre = '%s-%s.%s' % (default_exe,reffldr,IO.getFileName(fi))
        fiout = '%s/%s.%s' % (reffldr,outpre,fasta_ext)
        fionm = IO.getFileName(fiout)
        fi_nm = IO.getFileName(fi)        
        if os.path.isfile(fiout):
            log.write('Alignment (%s, %s) already done.' % (reffldr,fi_nm,rco))
            exe.assertDone(reffldr,fiout)
            log.incrementTimer()
            continue # already done
        cmd = 'echo `%s -i %s %s -o %s/%s -F fasta -r ON` > %s/%s.out' \
            % (exe.cmd,ref,fi,reffldr,outpre,reffldr,outpre)
        exe.add(cmd,fiout,fionm,reffldr,fi_nm)
    
    # Run commands.
    
    exe.run()