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
DEPNDS  = ['scipy','numpy']
LINKS   = ['http://github.com/AlexSafatli/LabBlouinTools/tarball/master#egg=package-1.0']

setup(name='abermusa',version=VERSION,description=DESCRIP,url=URL,author=AUTHOR,author_email=EMAIL,license='MIT',packages=['abermusa'],install_requires=DEPNDS,dependency_links=LINKS,zip_safe=False)