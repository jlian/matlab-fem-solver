# FEM Solver with MATLAB

Please read the PDF for more details.

This solver is not generic and solves just one problem. However, it should be easy to adapt it to different problems.

## How to use

Run the FreeFEM++ (examplemesh.edp) file in `/src` and to obtain the mesh file. It's also included pre-compiled in the same folder.

Then run the `femsolver.m` MATLAB file for analysis, and output.

To change the number of elements used for analysis, change the  `m` value in the FreeFEM++ file, and rerun.
