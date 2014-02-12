#!/usr/bin/env python

# timeGrabber.py
# -------------------------
# Jul 12, 2013; Alex Safatli
# -------------------------
# Grab the elapsed times from a
# set of /usr/bin/time output 
# for MStA runs and put
# them all into a single file.

import sys, os
from utils import IO, FASTAnet, PDBnet

times = {}
numst = {}

for fi in sys.argv[1:]:
    if os.path.isfile(fi) and fi != 'time.elapsed':
        # Get elapsed time.
        o = open(fi)
        read = o.read()
        frag = read.split('elapsed')[0]
        time = frag.split()[-1]
        o.close()
        times[fi] = time
        # Get number of structures.
        fipre = IO.getFileName(fi)
        fasta = fipre + '.fasta'
        logfi = fipre + '.log'
        pdbfi = fipre + '.pdb'
        if os.path.isfile(fasta):
            # Multiple structure alignment.
            f = FASTAnet.FASTAstructure(fasta)
            numst[fi] = len(f.orderedSequences)
        elif os.path.isfile(logfi): 
            # Damastes alignment?
            inp_cnt = 0
            o = open(logfi)
            for li in o:
                if 'Input: ' in li: inp_cnt += 1
            o.close()
            numst[fi] = inp_cnt          
        else:
            # Just set it to 0 for now.
            print 'NOTE; %s did not possess a .fasta or .log file.' % (fipre)
            numst[fi] = 0
    elif fi == 'time.elapsed': print 'Skipping %s; is output from this program.' % (fi)
    else: print 'Skipping %s; is not a file.' % (fi)

# Write to file.

o = open('time.elapsed','w')
o.write('File\tNum.Str.\tElapsed\n')
for time in times: o.write('%s\t%d\t%s\n' % (IO.getFileName(time),numst[time],times[time]))
o.close()
