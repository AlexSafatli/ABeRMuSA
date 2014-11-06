''' tmalign.py / Jul 25th, 2013; Alex Safatli
Plugin for ABeRMuSA adding support for the TM-Align pairwise aligner executable. '''

# Imports

import os
from labblouin import IO, FASTAnet, PDBnet

postprocess = True
default_exe = 'TMalign'
fasta_ext   = 'fasta'
ref_pos     = 1

def TMtoPDB(fi,out):
    '''
    Read a sup_atm output file from TMAlign and
    output to a proper PDB format.
    '''
    p = PDBnet.PDBstructure()
    o = open(fi)
    for li in o:
        if not li.startswith('ATOM') or len(li) < 47: continue
        resid = li[22:26].strip() # Residue ID
        chain = li[21:23].strip() # Chain
        if not chain in p.GetChainNames(): p.NewChain(chain)
        if not p.chains[chain].GetResidueByIndex(resid):
            resname = li[17:21].strip()
            res = PDBnet.PDBresidue(resid,resname)
            p.GetChain(chain).AddResidue(res)
        else: res = p.chains[chain][resid]
        serial = int(li[5:12].strip())
        atomn  = li[13:17].strip()
        xyz = li[26:].split()
        if len(xyz) != 3:
            x = float(li[31:38].strip())
            y = float(li[39:46].strip())
            z = float(li[47:].strip())
        else:
            x,y,z = xyz
            x = float(x)
            y = float(y)
            z = float(z)
        atom = PDBnet.PDBatom(serial,atomn,x,y,z,0,0,'','')
        res.AddAtom(atom)
    o.close()
    p.WriteFile(out)

def postProcess(fi,ref,log):
    
    reffldr = IO.getFileName(ref)
    finame  = IO.getFileName(fi)    
    
    # Parse output and acquire FASTA sequence.
    
    out_file = os.path.join(reffldr,finame + '.out')
    fas_file = os.path.join(reffldr,finame + '.' + fasta_ext)
    if not os.path.isfile(out_file): return
    o = open(out_file)
    l = o.read()
    o.close()
    find = l.rfind('residues)') # Signature denoting successful completion.
    if find == -1:
        log.write('tmalign: Alignment failed for <%s>.' % (finame))
        return
    raw = l.split('\n')
    lsindex = -2
    if raw[-1].strip() != '': lsindex = -1 # In case last line is not empty (blame TMalign).
    fasta = FASTAnet.FASTAstructure(uniqueOnly=False)
    fasta.addSequence('tmalign::target',raw[-4].strip())
    fasta.addSequence('tmalign::reference',raw[lsindex].strip())
    fasta.writeFile(fas_file)

    # Generate PDB file.
    
    pdfile  = os.path.join(reffldr,finame + '.sup_all_atm')
    if os.path.isfile(pdfile):
        TMtoPDB(pdfile,fi)
        log.write('tmalign: PDB file <%s> generated.' % (fi))
    else: log.write('tmalign: No Rasmol script <%s> found. Cannot generate PDB.' % (pdfile))

def executeCmd(args,ref,exe,log):
    
    # Decide on a reference folder.
    
    reffldr = IO.getFileName(ref)
    
    # Loop and add commands to exewrapper.
    
    for fi in [x for x in args if x != ref]:
        # For every file except the reference.
        outpre = '%s-%s.%s' % (default_exe,reffldr,IO.getFileName(fi))
        fiout = '%s/%s.%s' % (reffldr,outpre,fasta_ext)
        pdout = '%s/%s.%s' % (reffldr,outpre,'sup')
        fionm = IO.getFileName(fiout)
        fi_nm = IO.getFileName(fi)
        if os.path.isfile(fiout):
            log.writeTemporary('Alignment (%s, %s) is already done.' % (reffldr,fi_nm))
            exe.assertDone(reffldr,fiout)
            log.incrementTimer()
            continue # already done
        cmd = 'echo "`%s %s %s -o %s`" 2> %s > %s/%s.out' % (
            exe.cmd,fi,ref,pdout,fiout,reffldr,outpre)
        exe.add(cmd,fiout,fionm,reffldr,fi_nm)
    
    # Run commands.    
    
    exe.run()
        
    