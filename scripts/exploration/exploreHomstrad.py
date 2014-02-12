#!/usr/bin/env python

# exploreHomstrad.py; May 22, 2013
# usage: %prog dbloc
# Alex Safatli

import sys, os, glob, cPickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from utils import homstrad, homology
from utils.homology import writeSOAPSSAlignment

# Auxilary Functions

def hasNoFolders(path):
    for fi in glob.glob(os.path.join(path,'*')):
        if os.path.isdir(fi): return False
    return True

def getDirectories(path):
    fldrs = []
    for fi in glob.glob(os.path.join(path,'*')):
        if os.path.isdir(fi) and hasNoFolders(fi):
            h = homstrad.homstradFolder(fi)
            if len(h.getNames()) > 1:
                fldrs.append(h)
        elif os.path.isdir(fi):
            fldrs.extend(getDirectories(fi))
    return fldrs

# Script Body

if len(sys.argv) != 2:
    print 'usage: python homstradwalker.py homstraddbpath'
    exit()

dbloc = sys.argv[1]
pvals_soapsse = {}
pvals_matt    = {}
rawpv_soapsse = []
rawpv_matt    = []
rrmsds_homstr = []
align_lengths = []

# Get all Homstrad directories.
print '[....] Reading Homstrad database directories...'
fldrs = getDirectories(dbloc)

# Do all-for-all comparisons if not done already.
for fldr in fldrs:
    print '[%d/%d] Descending into %s...' % (fldrs.index(fldr)+1,len(fldrs),fldr.path)
    names = fldr.getNames()
    align_lengths.append(fldr.getAlignmentLength())
    if not os.path.isfile('homstrad.rrmsd'):
        alignedPDB = fldr.getAlignedPDB()
        alignedFAS = fldr.getFASTA()
        try:
            rr, rm = homology.rrmsd(alignedPDB,alignedFAS,rmsd=True)
            print '\trrmsd: Found RMSD %f and RRMSD %f...' % (rm,rr)
            rrmsds_homstr.append((rr,rm))
        except:
            print '\trrmsd: Could not find RMSD.'
            rrmsds_homstr.append((-1,-1))
    for name1 in names:
        for name2 in [x for x in names if x != name1]:
            if not os.path.isfile('homstrad.soapsse'):
                if name1 not in pvals_soapsse:
                    pvals_soapsse[name1] = {}                
                fast = fldr.writeFASTA('%s_%s.fasta' % (name1,name2),[name1,name2])
                alnf = '%s_%s.aln' % (name1,name2)
                pout = '%s_%s.pval' % (name1,name2)
                pdb1, pdb2 = [x for x in fldr.getPDBs() if name1 in x or name2 in x]
                try: writeSOAPSSAlignment(fast,pdb1,pdb2,alnf)
                except: 
                    os.remove(fast)
                    continue
                soapscmd = 'soapsse -a %s -b %s -n %s > %s' % (pdb1,pdb2,alnf,pout)
                print '\tsoapsse: Running %s vs. %s...' % (name1,name2)
                os.system(soapscmd)
                soaf = open(pout)
                pval = soaf.read().strip()
                soaf.close()
                try: pval = float(pval)
                except: pval = -1.0
                pvals_soapsse[name1][name2] = pval
                rawpv_soapsse.append(pval)
                print '\tsoapsse: Found P-value: %f...' % (pval)
                # Clean temporary files.
                os.remove(fast)
                os.remove(alnf)
                os.remove(pout)    
            if not os.path.isfile('homstrad.matt'):
                if name1 not in pvals_matt: pvals_matt[name1] = {}
                pdb1, pdb2 = [x for x in fldr.getPDBs() if name1 in x or name2 in x]
                mattcmd = 'Matt %s %s -o temp > temp.stdout' % (pdb1,pdb2)
                print '\tMatt: Running %s vs. %s...' % (name1,name2)
                os.system(mattcmd)
                p = -1.0
                try:
                    matf = open('temp.txt')
                    for line in matf:
                        if line.startswith('P-value'):
                            p = line.split()[-1].strip('\n')
                            break
                except: pass
                if p != -1.0:
                    try: p = float(p)
                    except: p = -1.0
                pvals_matt[name1][name2] = p
                rawpv_matt.append(p)
                print '\tMatt: Found P-value: %f...' % (p)
                os.system('rm temp.* > /dev/null 2>&1')

# Dump all p-value data to a pickle file and plot.
dump = 'homstrad.soapsse'
if len(pvals_soapsse) > 0:
    dumpf = open(dump,'w')
    cPickle.dump(pvals_soapsse,dumpf)
    dumpf.close()
    # Plot.
    plt.plot(range(len(rawpv_soapsse)),rawpv_soapsse)
    plt.title('Homstrad Database P-Value Distribution (SOAPSSe)')
    plt.ylabel('SOAPSSe P-Value')
    pp = PdfPages('homstrad_soapsse.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
else:
    dumpf = open(dump)
    pvals_soapsse = cPickle.load(dumpf)
    dumpf.close()
    for n in pvals_soapsse:
        for m in pvals_soapsse[n]:
            rawpv_soapsse.append(pvals_soapsse[n][m])
            
dump = 'homstrad.matt'
if len(pvals_matt) > 0:
    dumpf = open(dump,'w')
    cPickle.dump(pvals_matt,dumpf)
    dumpf.close()
    # Plot.
    plt.plot(range(len(rawpv_matt)),rawpv_matt)
    plt.title('Homstrad Database P-Value Distribution (Matt)')
    plt.ylabel('Matt P-Value')
    pp = PdfPages('homstrad_matt.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
else:
    dumpf = open(dump)
    pvals_matt = cPickle.load(dumpf)
    dumpf.close()
    for n in pvals_matt:
        for m in pvals_matt[n]:
            rawpv_matt.append(pvals_matt[n][m])
dump = 'homstrad.rrmsd'
if len(rrmsds_homstr) > 0:
    dumpf = open(dump,'w')
    cPickle.dump(rrmsds_homstr,dumpf)
    dumpf.close()
else:
    dumpf = open(dump)
    rrmsds_homstr = cPickle.load(dumpf)
    dumpf.close()

# If does not already exist, plot both against eachother.
if not os.path.isfile('homstrad.pdf'):
    plt.plot(rawpv_matt,rawpv_soapsse,'go')
    plt.title('Homstrad Database P-Value (Matt) vs. P-Value (SOAPSSe)')
    plt.ylabel('SOAPSSe P-Value')
    plt.xlabel('Matt P-Value')
    pp = PdfPages('homstrad.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
    
# If it does not already exist, plot RMSD vs RRMSD.
if not os.path.isfile('homstrad_rrmsd.pdf'):
    # Plot.
    rrvect, rmvect = [], []
    for it in rrmsds_homstr:
        rr, rm = it
        rrvect.append(rr)
        rmvect.append(rm)
    plt.plot(rrvect,rmvect,'go')
    plt.title('Homstrad Database RRMSD vs. RMSD Distribution')
    plt.xlabel('RRMSD')
    plt.ylabel('RMSD')
    pp = PdfPages('homstrad_rrmsd.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()    

# If it does not already exist, plot RRMSD vs lengths of alignment.
rrvect = []
for it in rrmsds_homstr:
    rr, _ = it
    rrvect.append(rr)
if not os.path.isfile('homstrad_rrlen.pdf'):
    # Plot.
    plt.plot(align_lengths, rrvect,'go')
    plt.title('Homstrad Database Alignment Lengths vs. RRMSDs')
    plt.ylabel('RRMSD')
    plt.xlabel('Alignment Length (# Residues)')
    pp = PdfPages('homstrad_rrlen.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
if not os.path.isfile('homstrad_rrlen_under200.pdf'):
    toplot = {'x':[],'y':[]}
    for r in xrange(len(rrvect)):
        if rrvect[r] >= 0 and align_lengths[r] < 200:
            toplot['x'].append(align_lengths[r])
            toplot['y'].append(rrvect[r])
    # Plot.
    plt.plot(toplot['x'],toplot['y'],'go')
    plt.title('Homstrad Database Alignment Lengths vs. RRMSDs if Length < 200')
    plt.ylabel('RRMSD')
    plt.xlabel('Alignment Length (# Residues)')
    pp = PdfPages('homstrad_rrlen_under200.pdf')
    plt.savefig(pp,format='pdf')
    pp.close()
