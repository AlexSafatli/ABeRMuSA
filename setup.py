''' An installation script for ABeRMuSA. '''

# Date:   Oct 17 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from setuptools import setup
from ABeRMuSA.ABeRMuSA import VERSION

DESCRIP = 'Multiple protein structural alignment leveraged by multiple pairwise alignments.'
URL     = 'http://www.github.com/AlexSafatli/ABeRMuSA'
AUTHOR  = 'Alex Safatli'
EMAIL   = 'safatli@cs.dal.ca'
DEPNDS  = []
LINKS   = []

setup(name='ABeRMuSA',version=VERSION,description=DESCRIP,url=URL,author=AUTHOR,author_email=EMAIL,license='MIT',packages=['abermusa'],install_requires=DEPNDS,dependency_links=LINKS,zip_safe=False)
