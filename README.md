ABeRMuSA
========

A Python application intended to perform multiple structural alignments on incredibly large numbers of protein structures that are not normally able to be done using traditional methods. Utilizes a strategy that leverages multiple pairwise alignments against a single reference.

ABeRMuSA stands for "Approximate Best Reference Multiple Structure Alignment" and alludes to the fact that a best reference is chosen through an exhaustive or approximate method and aligned in a pairwish fashion to all other structures in a given set of structures in order to come at an alignment.

Installation
--------

Running the following command while in this directory should install the ABeRMuSA executable to your /usr/bin/ directory.

  python setup.py install

Running
--------

  ABeRMuSA protein1.pdb protein2.pdb ... proteinN.pdb

A number of options are available, including initiating the heuristic (rather than exhaustive) search for a best reference by using the "-q" option and a number of iterations (or different references to try). You can always type "-h" in your command when executing the program in order to see all possible options.



