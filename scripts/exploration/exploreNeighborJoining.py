#!/usr/bin/env python

# exploreNeighborJoining.py
# -------------------------
# Jun 26, 2013; Alex Safatli
# -------------------------
# Explore the neighbor-joining
# landscape for a given set of
# reference folders (e.g., the
# full set of all references
# as found by an exhaustive search
# for a given dataset of PDBs
# by the Damastes aligner).
#
# Input: The Damastes XML file
# output for a run containing
# all files processed.
#
# usage: %prog XML_FILE

import sys, os, glob, cPickle, tempfile, igraph
from utils import PDBnet, logfile, IO
from Bio import Phylo
from itertools import combinations
from numpy import median, mean
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

for xmlfile in sys.argv[1:]:
    
    # Prepare temporary files.
    
    njfasta = tempfile.NamedTemporaryFile()
    njoutfi = tempfile.NamedTemporaryFile()
    njouttr = tempfile.NamedTemporaryFile()
    
    # Open XML file.
    
    filelist, bestref = [], None
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
        for ch in p.orderofchains: njfasta.write('>%s-%s\n%s\n' % (f,ch,p.ChainAsFASTA(ch)))
    njfasta.flush()
    os.fsync(njfasta)
        
    # Run Muscle.
    
    run = 'muscle -quiet -maxiters 1 -diags -sv -distance1 kbit20_3 -in %s -out %s -tree1 %s'\
        % (njfasta.name,njoutfi.name,njouttr.name)
    print 'MUSCLE; Running %s' % (run)
    os.system(run)
    
    # Get distance matrix.
    
    try: tree = Phylo.read(njouttr,'newick')
    except: continue
    dmat = {}
    for p1,p2 in combinations(tree.get_terminals(),2):
        d = tree.distance(p1,p2)
        if not p1.name in dmat: dmat[p1.name] = {}
        dmat[p1.name][p2.name] = d
    
    # Close files.
    
    njfasta.close()
    njoutfi.close()
    njouttr.close()
    
    # Get distance rank vector.
    
    dist_matrix = dmat
    name_matrix = {}
    norm_vector = {}
    rmsd_vector = {}
    rrmsd_vector = {}
    dist_means  = {}
    dist_trans  = {}
    for dli in dist_matrix:
        vals = [dist_matrix[dli][x] for x in dist_matrix[dli]]
        dist_means[dli] = mean(vals)
        norm_vector[dli[:-2]] = mean(vals)
        name_matrix[dli[:-2]] = {}
        for fli in dist_matrix[dli]: name_matrix[dli[:-2]][fli[:-2]] = dist_matrix[dli][fli]
    flat = [dist_means[x] for x in dist_means]
    mea, fmax = mean(flat), max(flat)
    for dli in norm_vector: norm_vector[dli] /= fmax
    
    # Go into every reference folder and determine
    # that reference's average RMSD, RRMSD.
    
    for dli in norm_vector:
        reffldr = IO.getFileName(dli)
        avgr, avgd, i = 0, 0, 0
        if not os.path.isdir(reffldr) or not os.path.isfile('%s/ref.pickl' % (reffldr)):
            continue
        for fi in glob.glob(os.path.join(reffldr,'*.pickl')):
            if 'ref.pick' in fi: continue
            o = open(fi)
            sc,rc,pv = cPickle.load(o)
            o.close()
            avgr += rc
            avgd += sc
            i += 1
        if (i == 0): continue
        rmsd_vector[dli] = avgr/i
        rrmsd_vector[dli] = avgd/i
        
    # Plot RMSD vs. Distances, RRMSD vs. Distances.
    
    name = IO.getFileName(xmlfile)
    x, y = [], []
    for dli in rmsd_vector:
        if bestref and dli == bestref: continue
        y.append(rmsd_vector[dli])
        x.append(norm_vector[dli])
    plt.plot(x,y,'go')
    if bestref and bestref in norm_vector:
        plt.plot([norm_vector[bestref]],[rmsd_vector[bestref]],'r*')
    flat = [norm_vector[j] for j in norm_vector]
    plt.axvline(x=mean(flat),color='r')
    plt.axvline(x=median(flat),color='b')
    plt.title('Average RMSD vs. Average Distance for %s' % (name))
    plt.xlabel('Average Genetic Distance')
    plt.ylabel('RMSD')
    pp = PdfPages('%s-rmsd.pdf' % (name))
    plt.savefig(pp,format='pdf')
    pp.close()
    
    x, y = [], []
    for dli in rmsd_vector:
        if bestref and dli == bestref: continue
        y.append(rrmsd_vector[dli])
        x.append(norm_vector[dli])
    plt.cla()
    plt.plot(x,y,'go')
    if bestref and bestref in norm_vector:
        plt.plot([norm_vector[bestref]],[rrmsd_vector[bestref]],'r*')
    plt.axvline(x=mean(flat),color='r')
    plt.axvline(x=median(flat),color='b')
    plt.title('Average RRMSD vs. Average Distance for %s' % (name))
    plt.xlabel('Average Genetic Distance')
    plt.ylabel('RRMSD')
    pp = PdfPages('%s-rrmsd.pdf' % (name))
    plt.savefig(pp,format='pdf')
    pp.close()
    
    print 'DONE; PDFs of plots found in current directory.'