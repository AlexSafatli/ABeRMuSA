''' mustang.py; May 10th, 2013; Alex Safatli
Plugin for ABeRMuSA adding support for the MUSTANG pairwise aligner executable. '''

# Imports

import os
from labblouin import IO

postprocess = False
default_exe = 'mustang'
fasta_ext   = 'afasta'
ref_pos     = 0

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
            log.writeTemporary('Alignment (%s, %s) already done.' % (reffldr,fi_nm,rco))
            exe.assertDone(reffldr,fiout)
            log.incrementTimer()
            continue # already done
        cmd = 'echo `%s -i %s %s -o %s/%s -F fasta -r ON` > %s/%s.out' \
            % (exe.cmd,ref,fi,reffldr,outpre,reffldr,outpre)
        exe.add(cmd,fiout,fionm,reffldr,fi_nm)
    
    # Run commands.
    
    exe.run()