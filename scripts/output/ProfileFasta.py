# ProfileFasta.py
# ----------------------------------
# Add a new sequence to an alignment.

import sys,os
from glob import glob as g
from subprocess import Popen
from utils.FASTAnet import FASTAstructure as F

def  reformat_fasta(fastafile,multiple=False):
    f=open(fastafile)
    s=''
    for line in f:
        if line.startswith('>'):
            if multiple:
            
                name=line.replace('>','').strip()
                name = name[:name.find(':')+1]
                pdb=name[:4]
                nl='>'+name+pdb+'\n'
            else:
                nl = line
            s+=nl
        else:
            s+=line
    fout=open(fastafile,'w')
    fout.write(s)
    fout.close()

if __name__ == '__main__':

    # If not imported.

    prefix=sys.argv[1]
    if '-m' in sys.argv[2:]: multiple=True
    else: multiple=False
        
    files=g('*.fasta')
    oal=Popen('muscle -profile -in1 %s -in2 %s -out currental.afa -quiet -verbose'
              %(files[0],files[1]),shell=True)
    oal.wait()
    for f in files[2:]:
        print 'Including',f,'in the overall alignment.'
        nal=Popen('muscle -profile -in1 currental.afa -in2 %s -out currental.afa -quiet -verbose'
                  %(f),shell=True)
        nal.wait()
    
    os.rename('currental.afa',prefix+'.fasta')
    
    ali=F(prefix+'.fasta')
    print 'Removing duplicate and redundant sequences.'
    ali.writeFile(prefix+'.fasta')

