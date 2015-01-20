ABeRMuSA
========

A Python application that is capable of performing multiple structural alignments on incredibly large numbers of three-dimensional protein structures that are not normally able to be done using traditional methods. Utilizes a strategy that leverages multiple pairwise alignments against a single reference.

ABeRMuSA stands for "**A**pproximate **Be**st **R**eference **Mu**ltiple **S**tructure **A**lignment", alluding to the fact that a best reference is chosen through an exhaustive or approximate method and aligned in a pairwish fashion to all other structures in a given set of structures in order to develop an alignment.

Installation
--------

Running the following command while in this directory should install the ABeRMuSA executable to your /usr/bin/ directory.

    python setup.py install

If you have any problems installing the software, [e-mail the author](mailto:safatli@cs.dal.ca) or submit an issue.

Running
--------

    ABeRMuSA protein1.pdb protein2.pdb ... proteinN.pdb

A number of options are available, including initiating the heuristic (rather than exhaustive) search for a best reference by using the "-q" option and a number of iterations (or different references to try). You can always type "-h" in your command when executing the program in order to see all possible options. Typical behavior of the software will output a final PDB alignment file, in addition to a number of other auxiliary files such as a FASTA alignment and logfile, with a prefix of the user's choice.

