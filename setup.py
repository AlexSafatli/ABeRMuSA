''' An installation script for ABeRMuSA. '''

# Date:   Oct 17 2014
# Author: Alex Safatli
# E-mail: safatli@cs.dal.ca

from setuptools import setup
from abermusa.__version__ import VERSION

DESCRIP = 'Multiple protein structural alignment leveraged by multiple pairwise alignments.'
URL     = 'http://www.github.com/AlexSafatli/ABeRMuSA'
AUTHOR  = 'Alex Safatli'
EMAIL   = 'safatli@cs.dal.ca'
DEPNDS  = ['scipy','numpy']
SCRIPTS = ['ABeRMuSA']
LINKS   = ['http://github.com/AlexSafatli/LabBlouinTools/tarball/master#egg=package-1.0']
PKGDATA = {'abermusa':['plugins/*']}

setup(name='abermusa',version=VERSION,description=DESCRIP,url=URL,author=AUTHOR,author_email=EMAIL,license='MIT',packages=['abermusa'],package_data=PKGDATA,install_requires=DEPNDS,dependency_links=LINKS,scripts=SCRIPTS,zip_safe=False)
