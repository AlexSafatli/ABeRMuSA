#!/usr/bin/env python

# rerunFailed.py
# -------------------------
# Jul 22, 2013; Alex Safatli
# -------------------------
# Simple utility to go over
# a failed file output from
# summarizer.py and rerun Damastes
# on them based on a given shell
# script.

import sys, os

if len(sys.argv) == 1:
    print 'usage: %s shell_script failed_file output_file' % (sys.argv[0])
    exit()

shell_script = sys.argv[1]
failed_file  = sys.argv[2]
output_file  = sys.argv[3]
o = open(failed_file)
l = o.readlines()
o.close()
if os.path.isfile(output_file): os.remove(output_file)
for li in l: os.system('cat %s | grep %s/ >> %s' % (shell_script,li.strip(),output_file))