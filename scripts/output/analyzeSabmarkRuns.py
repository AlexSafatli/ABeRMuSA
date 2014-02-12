#!/usr/bin/env python

# analyzeSabmarkRuns.py; Jun 07, 2013
# usage: %prog sabdbfile reffldr [reffldr2] [reffldr3] ...
# Alex Safatli
# -------------------------------------------------------
# Read the reference folder output from a Sabmark run.
# -------------------------------------------------------

import sys, os, cPickle, numpy
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from scipy.stats import norm, normaltest
from utils import PDBnet

# Script Body /////////////////////////////////////////////////////////////////////

if len(sys.argv) < 3:
    print 'usage: %s sabdbfile reffldr [reffldr2] [reffldr3] ...' % (sys.argv[0])
    exit()

dict_rrmsd    = {}
raw_rrmsd     = []
plot_rrmsd    = []
plot_matt     = []
found_greater = []

# Assumes Matt.
for reffldr in sys.argv[2:]:
    # For every reference folder given from Damastes...
    pickl_file = os.path.join(reffldr,'ref.pickl')
    if not os.path.isfile(pickl_file):
        print '%s is not completely scored or not a ABeRMuSA reference folder.' % (reffldr)
        continue
    o = open(pickl_file) # Open pickle file.
    scores, _, _, _ = cPickle.load(o)
    matt_dict = {}
    rrmsd_dict = {}
    o.close()
    for s in scores:
        val = None
        o = open('%s/%s.txt' % (reffldr,s))
        for line in o:
            if line.startswith('P-value'):
                l = line.split(':')
                val = float(l[-1].strip())
                break
        o.close()
        o = open('%s/%s.pickl' % (reffldr,s))
        sc,rc,pv = cPickle.load(o)
        o.close()
        if val:
            if (val > 1):
                found_greater.append(('MATT ','%s/%s'%(reffldr,s),'pval: %f' % (val)))
                val = 1.0
            if (pv > 1):
                found_greater.append(('RRMSD','%s/%s'%(reffldr,s),'rrmsd: %s'%(sc),'pval: %f' % (pv)))
                pv = 1.0
            matt_dict[s] = val
            rrmsd_dict[s] = pv
            raw_rrmsd.append(sc)
            plot_rrmsd.append(rrmsd_dict[s])
            plot_matt.append(matt_dict[s])
        namesplit = s.split('.')
        ref, oth = namesplit[0].strip('Matt-'), namesplit[1]
        if not ref in dict_rrmsd: dict_rrmsd[ref] = {}
        dict_rrmsd[ref][oth] = (sc,rc,pv)

# Open Sabmark DB file and determine true and false positive nature
# of all files.

print '\nOpening database in order to determine truth values...\n'
true_false = {}
if not os.path.isfile(sys.argv[1]):
    print '%s not found.' % (sys.argv[1])
    exit()
fh = open(sys.argv[1])
for line in fh:
    l = line.split('\t')
    if len(PDBnet.PDBstructure(l[0].strip())) < 100: 
        print '%s < 100 residues; ignoring.' % (l[0])
        continue
    temp = l[0].split('.')[0]
    n = temp.split('/')[-1]
    true_false[n] = (l[1].strip() == 'True')

# Plot the histogram.
if len(found_greater) > 0:
    print '\nThe following P-values were found to be greater than 1; details are provided on each line.\n'
    for g in found_greater: print g
plt.hist([plot_matt,plot_rrmsd],color=['crimson','burlywood'],label=['Matt','RRMSD'])
plt.title('Matt P-Values vs. RRMSD P-Values by Frequency')
pp = PdfPages('hist.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.cla()

# Plot the P-values against eachother.
plt.plot(plot_matt,plot_rrmsd,'ro')
plt.title('Matt P-Values vs. RRMSD P-Values by Value')
pp = PdfPages('pvals.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.cla()

# Magic: find only false positives vs. true positives.
plot_tf = []
plot_tp = []
for r in dict_rrmsd:
    for o in dict_rrmsd[r]:
        if (not r in true_false) or (not o in true_false):
            print '<%s, %s> not in truth table; ignoring.' % (r,o)
            continue
        if (true_false[r] and true_false[o]):
            plot_tp.append(dict_rrmsd[r][o][0])
        if (not true_false[r] and true_false[o]) or \
           (true_false[r] and not true_false[o]):
            plot_tf.append(dict_rrmsd[r][o][0])

# Plot them.
plt.hist(plot_tf,bins=100)
y,binEdges=numpy.histogram(plot_tp,bins=100)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.plot(bincenters,y,'r-')
plt.hist(plot_tp,bins=100,color='red')
_, pv = normaltest(plot_tf)
mean, stdv = norm.fit(plot_tf)
print '\nNormality Test (true vs. false): mean: %f, stdv: %f, pval: %f' % (mean,stdv,pv)
_, pv = normaltest(plot_tp)
mean, stdv = norm.fit(plot_tp)
print 'Normality Test (only true positive): mean: %f, stdv: %f, pval: %f' % (mean,stdv,pv)
plt.title('RRMSD Values for True vs. False Positive Alignments')
pp = PdfPages('tfrrmsd.pdf')
plt.savefig(pp,format='pdf')
pp.close()
plt.cla()


