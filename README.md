# parallel-Goldstein
Parallel Goldstein phase unwrapping algorithm based on breadth first search

This repository contains the source code and data files of the parallel golstein algorithm for phase unwrapping published in https://doi.org/10.1016/j.optlaseng.2019.105827

See the LICENSE file for the licensing statement and Disclaimer of Warranty.

The source code was written in C on gcc compiler using Linux Ubuntu operating system. The organization of the directory qgpu is the following:

1. lib, include, src: Implementation of the quality guide phase unwrapping algorithm with the pruning strategy.

2. data: Contains the wrapped phase maps, quality maps and mask files used in the submitted paper.
3. matlab: scripts for generating simulated phase distributions of any size

Compiling:

The implementation of the parallel Goldstein algorithm contains three directories, say include, lib and src. The header files in each implementation are in the include directory and the source files are in the src directory. To generate the executable program, locate in the corresponding src directory and type make in order to compile the project. To execute the program, type ./main in the current src directory. This implemetation requires to install the OpenMP library on your system.

The project can be compiled either in MAC OS x and in Linux operating systems with minor modifications.
