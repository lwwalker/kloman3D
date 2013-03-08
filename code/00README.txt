Updated 3/8/13

kloman3D: kernel learning of monoclonal antibody neutralization in 3D

This software is designed to construct a predictive kernel model for 
antibody neutralization taking into account the 3D structure of
the HIV env protein.  

Requires:
Python 2.x (developed on 2.7)
libsvm

Usage: kloman3D.py [optional flags] responsefile outfile
Manditory parameters:
responsefile -> filename of file containing training classes, should be 2 column comma-delimited file with column 0 being sample names matching fasta file, and column 1 being 1 or -1 indicating classes
outfile -> file to write output to

Optional parameters:
-M --matrix -> filename of pre-computed array of kernel matrices
-m --map -> filename of seperate alignment map of sequences to consensus
-f --fasta -> filename for fasta file with reference sequence (if passing map) or all sequences if no map - mandatory of inot loading pre-computed matrixes
-G --gp120 -> pdb filename for gp120 - mandatory if not loading pre-computed matrixes
-g --gp41 -> optional pdb filename for gp41
-N --norm -> normalize the kernel matrixes
-r --radius -> The "bubble" radius for the kernel (default 10 Angstroms), may pass a single value, or comma-delimited list
-s --stdev -> Standard deviation for the word kernel function (default 10 Angstroms)
-x --nocross -> use non-cross-validation SVMs for feature selection
-a --aminoacid -> filename for an amino acid partial-match matrix
-c --cSVM -> perform c-svm with the passed regularization parameter
-n --nuSVM -> perform nu-svm with the passed regularization parameter (default with nu = 0.2)
-v --saveKM -> Save the list of position-wise kernel matrixes to supplied filename for future use
-h --help -> Display this message

Examples:

./kloman3Dv05.py -m ../rawData/2F5.py -G ../supportFiles/3DNN.pdb -f ../supportFiles/HXB2_nogaps.fasta ../analysis/2F5.resp 2F5.out

./kloman3Dv05.py -M test.km ../analysis/2F5.resp test.out

