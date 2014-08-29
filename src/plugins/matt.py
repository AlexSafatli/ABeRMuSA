''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++

matt.py
May 6th, 2013; Alex Safatli
Plugin for ABeRMuSA adding support for the MATT pairwise aligner executable. '''

# Imports

import os
from utils import IO

postprocess = False
default_exe = 'Matt'
fasta_ext   = 'fasta'
ref_pos     = 1

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
            log.writeTemporary('Alignment (%s, %s) already done.' % (reffldr,fi_nm))
            exe.assertDone(reffldr,fiout)
            log.incrementTimer()
            continue # already done
        cmd = 'echo `%s %s %s -o %s/%s` > %s/%s.out' % (exe.cmd,fi,ref,reffldr,outpre,reffldr,outpre)
        exe.add(cmd,fiout,fionm,reffldr,fi_nm)
    
    # Run commands.    
    
    exe.run()
        
    