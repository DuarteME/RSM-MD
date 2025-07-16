# RSM-MD
 This is a Python script to calculate the reciprocal space map of Ga2O3 simulation cells obtained by molecular dynamics.

 Before running RSM_calculate.py, edit the name of the source files (in this case, they are of the form "{i}-frame.data") and the output files ("RSM_{i}-frame.dat"). The number of rows (i.e., particles) should also be given.

 Then, the desired reciprocal space region and step should be given in reciprocal Angstrom units. In this case, the discretisation of the reciprocal space is performed in 3 directions ({a*, b, c}), between -82 and 82 nm-1, with 51 points per direction. 

 In order to keep up with the usual conventions, the input is given in Angstrom (as used by LAMMPS) and the output in given in nm-1. 

 This is easily extended to other materials by using the appropriate Cromer-Mann coefficients.

 The results can be easily plotted by the companinon script RSM_plot.py.
