# exploreProcrustes.py
# ---------------------------------------------
# Jun 28, 2013; Alex Safatli
# ---------------------------------------------
# Test the viability of using Procrustes
# in the context of Damastes alignments.
#
# Input: The Damastes XML file
# output for a run containing
# all files processed.
#
# usage: %prog PROC_BIN XML_FILE

import sys, os, glob, tempfile, igraph, math
from utils import PDBnet, FASTAnet, logfile, IO
from Bio import Phylo
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# //////////////// FROM GMtoRMSD ////////////////

def GMtoMatrix(prefix):
    f = open(prefix + '.gm')
    s = {} ## s for structures
    for line in f:
        if line.strip().endswith(';'):
            line = line.strip()[:-1]
        temp = [[], [], []]
        bline = line[line.find(';')+1:].strip().split(';')
        for i in range(len(bline)):
            temp[i%3].append(float(bline[i]))
        s[line[:line.find(';')]] = temp
        # Check for error (in case parsing in this function is bad...)
        if not (len(temp[0]) == len(temp[1]) == len(temp[2])):
            print 'Parsing error for structure ' + line[:line.find(';')]
    return s


def RMSDfromMEAN(p,m):
    # This code was taken from StructureOutlier because that file couldn't be imported!!!

    # Calculates the RMSD for the meanshape (m) and the problem shape (p)
    sum_dist=0.0
    count = 0
    for j in range(len(m[0])):
        d =(m[0][j]-p[0][j])**2 + (m[1][j]-p[1][j])**2 + (m[2][j]-p[2][j])**2
        sum_dist += d
        count += 1
    # calculate the sqrt of the mean deviation
    RMSD = math.sqrt(sum_dist/count)

    return RMSD

# /////////////// SCRIPT BODY ///////////////////

bins = {}
results = {}

for xmlfile in sys.argv[2:]:

    # Prepare temporary files.

    njfasta = tempfile.NamedTemporaryFile()
    njoutfi = tempfile.NamedTemporaryFile()

    # Open XML file.

    filelist, pdblist, bestref, bestind = [], {}, None, -1
    xml = logfile.XMLfile(xmlfile)
    xml.read()
    for child in xml.root:
        if child.tag == 'file':  filelist.append(child.text)
        elif child.tag == 'reference': bestref = child.text
    if not bestref: print 'WARNING; Best reference was not found in XML file given.'
    else: print 'REFERENCE; %s' % (bestref)

    # Extract FASTA files from PDBs.

    for f in filelist:
        # Will be a PDB file. Add it to input fasta.
        p = PDBnet.PDBstructure(f)
        pdblist[f] = p
        for ch in p.orderofchains: njfasta.write('>%s\n%s\n' % (f,p.ChainAsFASTA(ch)))
    njfasta.flush()
    os.fsync(njfasta)

    # Run Muscle.

    run = 'muscle -quiet -maxiters 1 -diags -sv -distance1 kbit20_3 -in %s -out %s'\
        % (njfasta.name,njoutfi.name)
    print 'MUSCLE; Running %s' % (run)
    os.system(run)

    # Get output alignment.

    f = FASTAnet.FASTAstructure(njoutfi.name)
    seqs = f.sequences
    cons = range(len(f.orderedSequences[0].sequence))
    for seq in seqs:
        seq = seqs[seq].sequence
        for i in range(len(seq)):
            if seq[i] == '-' and i in cons: cons.remove(i)

    # For conserved regions, get alpha carbons from PDB files.

    prefi = '%s_sequence' % (IO.getFileName(xmlfile))
    outgm = '%s.gm' % (prefi)
    print 'GMFILE; Writing to %s.' % (outgm)
    gmfile = open(outgm,'w')
    i = 1
    for p in pdblist:
        if not p in seqs:
            print 'WARNING; %s not found in output from Muscle.' % (p)
            continue        
        gmfile.write('>%s;' % (p))
        print 'GMFILE; %s ("%d")' % (p,i)
        pdb = pdblist[p]
        ind = 0
        ch = pdb.orderofchains[0]
        inds = pdb.chainsOrder[ch]
        for it in range(len(seqs[p].sequence)):
            if seqs[p].sequence[it] != '-' and it in cons:
                ca = pdb.chains[ch][inds[ind]].GetCA()
                if ca: gmfile.write('%f;%f;%f;' % (ca.x,ca.y,ca.z))
                else: print 'WARNING; No alpha carbon at %d (%s).' % (it,p)
                ind += 1
            elif seqs[p].sequence[it] != '-': ind += 1
        gmfile.write('\n')
        if bestref and p == bestref: bestind = i
        i += 1
    gmfile.close()

    # Run the Procrustes analysis.

    proc_bin = sys.argv[1]
    proc_cmd = 'python %s %s 3 -scale' % (proc_bin,prefi)
    print 'PROCRUSTES; Running %s' % (proc_cmd)
    os.system(proc_cmd)

    # Compare to mean shape and get vector.

    vector = {}
    try:
        meanshape = GMtoMatrix('mshape')
        meanshape = meanshape[meanshape.keys()[0]]
        full      = GMtoMatrix(prefi)
    except:
        print 'ERROR; Could not proceed with RMSD calculations.'
        njfasta.close()
        njoutfi.close()
        continue
    names = full.keys()
    for name in names: vector[name] = RMSDfromMEAN(full[name],meanshape)
    
    # Sort vector and return sorted list of names.
    
    vecfile = '%s.seq_vector' % (IO.getFileName(xmlfile))
    vecfh = open(vecfile,'w')    
    keys = sorted(vector,key=lambda d:vector[d])
    i = 0
    for key in keys:
        out = 'STRUCTURE; %s (index: %d, rmsd: %f)' % (key,keys.index(key),vector[key])
        print out
        vecfh.write(out + '\n')
        ind = int(key.strip('"'))
        if bestref and ind == bestind:
            if not i in bins: bins[i] = 0
            bins[i] += 1
            results[xmlfile] = i
        i += 1

    # Close files.

    vecfh.close()
    njfasta.close()
    njoutfi.close()
    
# Plot counts of best reference.

bins_out = []
for x in bins:
    if x <= 2:
        for i in range(bins[x]): bins_out.append(x)
    else:
        for i in range(bins[x]): bins_out.append(3)
plt.hist(bins_out,bins=4)
pp = PdfPages('procr.pdf')
plt.savefig(pp,format='pdf')
pp.close()
for result in results:
    print 'ALIGNMENT (%s); best reference was at %dth position.' % (result,results[result])
print 'DONE.'