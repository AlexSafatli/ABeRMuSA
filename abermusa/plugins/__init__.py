''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++ '''

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

# DO NOT EDIT BELOW ++

import os

thisdir = os.path.split(os.path.realpath(__file__))[0]
itlist = os.listdir(thisdir)
__all__ = [os.path.split(x)[-1].strip('.py') for x in itlist if x.endswith('.py') and not x.endswith('__init__.py')]

# DO NOT EDIT ABOVE ++