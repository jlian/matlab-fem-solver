# FEM Solver with MATLAB

<img src="https://raw.githubusercontent.com/jlian/fem-solver/master/doc/matlab.png" width="400">

Please read the [PDF](https://github.com/jlian/fem-solver/blob/master/fem_solver.pdf) for more details.

This solver is not generic and solves just one problem. However, it should be easy to adapt it to different problems. Not tested with Octave.

## How to use

Run the FreeFEM++ (`examplemesh.edp`) file in `/src` and to obtain the mesh file. It's also included pre-compiled in the same folder.

Then run the `femsolver.m` MATLAB file for analysis, and output.

To change the number of elements used for analysis, change the  `m` value in the FreeFEM++ file, and rerun.
