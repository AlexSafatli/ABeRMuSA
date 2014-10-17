''' A plugin interface for ABeRMuSA. '''

# Plugins require the following function prototypes implemented:
# def executeCmd(prefix,args,ref,exe,log) --> executes command
#
# Plugins also require these variable definitions:
# default_exe = 'ExecutableName'
# fasta_ext   = 'fasta'

# DO NOT EDIT BELOW ++

import os

thisdir = os.path.split(os.path.realpath(__file__))[0]
itlist = os.listdir(thisdir)
__all__ = [os.path.split(x)[-1].strip('.py') for x in itlist if x.endswith('.py') and not x.endswith('__init__.py')]

# DO NOT EDIT ABOVE ++