# Minimal Solvers for Indoor UAV Positioning
This repository contains implementations from [1]. For more information, see the arXiv paper.

## Depedencies
The implementation uses Eigen 3 (older versions are not compatible), which is
a C++ template library for linear algebra: matrices, vectors,
numerical solvers, and related algorithms.

Installation for Ubuntu/Debian:

    $ apt-get install libeigen3-dev

Tested on version 3.3.4.

## Using the solvers in MATLAB
It is possible to MEX-compile the solver and use it in MATLAB.
In order to do so you may use the file `matlab/compile_mex.m`, which
should result in an executable.

Note that your local Eigen path may be different, e.g. `/usr/local/include/eigen3`.

Tested on version R2019a Linux (64-bit).

## Using the solvers in Python
This is unfortunately not supported yet; however, it is possible using eigency.

## Citations
Please cite us if you use the code:
[1] Valtonen Örnhag et al. (2020), Minimal Solvers for Indoor UAV Positioning.
```
@misc{valtonenoernhag-etal-2020-arxiv,
    title={Minimal Solvers for Indoor UAV Positioning},
    author={Marcus Valtonen~{\"O}rnhag and Patrik Persson and M{\aa}rten Wadenb{\"a}ck and Kalle {\AA}str{\"o}m and Anders Heyden},
    year={2020},
    eprint={--},
    archivePrefix={arXiv},
    primaryClass={cs.CV}
}
```
