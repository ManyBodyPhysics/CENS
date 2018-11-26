These folders contain (parallel versions will made available later) codes to renormalize a nucleon-nucleon interaction as well as computing an effective interaction for effective Hilbert spaces typically
used in so-called shell-model calculations in nuclear physics. See the tutorial folder for more information.

If you move into the folder  Benchmarks/VeffExample/Calcium40  you will find 
benchmarks for the full renormalized force and an effective interaction calculation
for the fp shell.
The input file for ca40 is set up for  6 major shells.
The file bhf.ini  defines the inputs to the effective interaction calculations.
The output files for the effective interactions are tailored to the fp-shell only
using the N3LO interaction and many-body perturbation theory to second order. The matrix elements for the N3LO interaction are included here.
The interaction matrix elements have been generated with the renomalization codes using the inputs defined in the 
renorm.ini file.
The effective interaction files are generated with the bhf.ini file. No executable are listed in this folder. These would have to be generated again using the veffective codes and the vrenormalization codes available from this website.
The codes were compiled using gfortran version 5.3.0 (Homebrew gcc 5.3.0)  on a Macbook PRO 3 GHz Intel Core i7
on May 25 2016
