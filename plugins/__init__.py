# This folder is a "package". All python files inside are considered
# to be a part of this package. To import this directory, use the
# following syntax in a parent directory:
#
# from FOLDERNAME import *
#
# Plugins require the following function prototypes implemented:
# def executeCmd(prefix,args,ref,exe,log) --> executes command
#
# Plugins also require these variable definitions:
# default_exe = 'ExecutableName'
# fasta_ext   = 'fasta'

import os

thisdir = os.path.split(os.path.realpath(__file__))[0]
itlist = os.listdir(thisdir)
__all__ = [os.path.split(x)[-1].strip('.py') for x in itlist if x.endswith('.py') and not x.endswith('__init__.py')]
