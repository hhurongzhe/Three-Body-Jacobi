# Three-Body-Jacobi

Calculating matrix elements of chiral 3N interactions under LS-coupled three-body partial-wave basis (in Jacobi coordinates), written in C++.

### Requirements

- This code requires openmp only.
- A C++ compiler with a makefile would be enough. So I provide an exmaple in “Makefile”.

### Usage

- use Mathematica notebooks to generate codes in ./deps-mma
- copy the generated .out files into src/aPWD3_part_*.cpp
- set N_channels in src/main_part_*.cpp
- compile and run
- save channel information and (2J,P,2T) in a txt file

### Method

To do partial-wave projection for 3N interactions, I follow the "aPWD" method developed by J. Golak. See "A new way to perform partial wave decompositions of few-nucleon forces" if interested.

### Structure of Code

- src/*: codes.
- tool/\*.py: some simple python tools to generate mesh points, channels index and plot.
- xmake.lua: Xmake’s project description file.
- Makefile: you may use it to compile.
- data/: directory of storing 3BMEs.

### Remarks

This code is mainly designed for generating non-locally regulated 3N interaction in the Jacobi coordinate (up to N2LO). However, it is possible to extend for local 3N interaction or for N3LO terms.

If you have any needs or questions, just contact me: rongzhe_hu@pku.edu.cn or rongzhehuu@gmail.com !
