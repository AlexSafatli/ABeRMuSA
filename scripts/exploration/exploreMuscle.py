#!/usr/bin/env python

# exploreMuscle.py
# -------------------------
# Jun 18, 2013; Alex Safatli
# -------------------------
# Perform neighbor-joining on
# input sequences in addition
# to sequence alignment by Muscle.
#
# usage: %prog SH_FILE

import sys, os
from utils import PDBnet
from Bio import Phylo
from itertools import combinations
from numpy import median, mean
from matplotlib import pyplot as plt

# From quickref.py
#def distMatrix(self):
    #'''
    #Perform a Muscle sequence alignment on extracted
    #PDB sequences in order to retrieve a distance
    #matrix.
    #'''
    #njfasta = tempfile.NamedTemporaryFile()
    #njoutfi = tempfile.NamedTemporaryFile()
    #njouttr = tempfile.NamedTemporaryFile()
    #for f in self.filelist:
        ## Will be a PDB file; add it to input fasta.
        #p = PDBnet.PDBstructure(f)
        #for ch in p.orderofchains: njfasta.write('>%s-%s\n%s\n' % (f,ch,p.ChainAsFASTA(ch)))
    #njfasta.flush()
    #os.fsync(njfasta)        
    #run = 'muscle -quiet -maxiters 1 -diags -sv -distance1 kbit20_3 -in %s -out %s -tree1 %s'\
        #% (njfasta.name,njoutfi.name,njouttr.name)
    #os.system(run)
    ## Get distance matrix from file.
    #tree = Phylo.read(njouttr,'newick')
    #dmat = {}
    #for p1,p2 in combinations(tree.get_terminals(),2):
        #d = tree.distance(p1,p2)
        #if not p1.name in dmat: dmat[p1.name] = {}
        #dmat[p1.name][p2.name] = d
    #njfasta.close()
    #njoutfi.close()
    #njouttr.close()
    #return dmat
#def getDistRankVector(self):
    #'''
    #With a given distance matrix, rank all files
    #by their distance from the mean distance.
    #'''
    #dist_matrix = self.distMatrix()
    #dist_means  = {}
    #dist_trans  = {}
    #for dli in dist_matrix:
        #vals = [dist_matrix[dli][x] for x in dist_matrix[dli]]
        #dist_means[dli] = mean(vals)
    #mea = mean([dist_means[x] for x in dist_means])
    #for dli in dist_means: dist_trans[dli] = abs(dist_means[dli]-mea)
    #temp_vector = sorted(dist_trans,key=lambda x:dist_trans[x])
    #dist_vector = []
    #for t in temp_vector: dist_vector.append(t[:-2])
    #return dist_vector


if len(sys.argv) != 2: exit()

o = open(sys.argv[1])
r = o.read().split('\n')
o.close()

distmats = {}
actuals  = []
statact  = []

for line in r:

    # Read line.
    
    files = []
    if line.startswith('echo'): continue
    l = line.split()
    nexti = False
    nam = ''
    for i in l:
        if nexti:
            nexti = False
            nam = i
        elif i == '-p': nexti = True
        elif i.startswith('-') or 'ABeRMuSA' in i: continue
        files.append(i)

    infasta = 'nj-in.fasta' 
    
    # Look for best reference in log file.
    
    logf = nam + '.log'
    if not os.path.isfile(logf): continue
    o = open(logf)
    k = o.read().split('\n')
    o.close()
    found = None
    for i in k:
        if 'Best reference found to be' in i:
            i = i.split('(')
            found = i[1][:-2]
            break
    if not found: continue
    ref = found
    
    # Extract FASTA files from PDBs.
    
    for f in files:
        if os.path.isfile(f):
            # Will be a PDB file. Add it to input fasta.
            p = PDBnet.PDBstructure(f)
            o = open(infasta,'a')
            for ch in p.orderofchains:
                o.write('>%s-%s\n%s\n' % (f,ch,p.ChainAsFASTA(ch)))
            o.close()
        
    # Run Muscle.
    
    run = 'muscle -quiet -maxiters 1 -diags -sv -distance1 kbit20_3 -in %s -out nj-out.fasta -tree1 nj-out.phy' % (infasta)
    os.system(run)
    
    # Get distance matrix.
    
    th   = open('nj-out.phy')
    tree = Phylo.read(th,'newick')
    dmat = {}
    for p1,p2 in combinations(tree.get_terminals(),2):
        d = tree.distance(p1,p2)
        if not p1.name in dmat: dmat[p1.name] = {}
        dmat[p1.name][p2.name] = d
    th.close()
    
    # Print results to file.
    
    oh  = open('nj-out.dist','w')
    for d in dmat: oh.write('%s\t%s\n' % (d,dmat[d]))
    oh.close()
    
    # Report stats.
    
    dli  = [dmat[x] for x in dmat]
    means = []
    flat = []
    for d in dli:
        li = [d[x] for x in d]
        means.append(mean(li))
        flat.extend(li)
    ma = max(flat)
    i = min(flat)
    e = mean(means)
    print 'UNNORMALIZED: Maximum Distance: %f, Minimum Distance: %f, Mean Distance: %f' % (ma,i,e)
    
    # Normalize.
    norm = [x/ma for x in flat]
    m = max(norm)
    i = min(norm)
    e /= ma
    print 'NORMALIZED  : Maximum Distance: %f, Minimum Distance: %f, Mean Distance: %f' % (m,i,e)
    distmats[nam] = dmat
    
    # Find reference structure.
    
    found = None
    for name in dmat:
        if ref in name:
            found = dmat[name]
            break
    if found:
        mea = mean([found[x] for x in found])
        print '%s found at distance: %f' % (ref,mea/ma)
        actuals.append(mea/ma)
        statact.append(e)

# Plot.

plt.plot(actuals,statact)
plt.show()