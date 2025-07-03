# Three-Body-Jacobi

Calculating matrix elements of chiral 3N interactions under LS-coupled three-body partial-wave basis (in Jacobi coordinates), written in C++17.

### Requirements

- This code requires openmp only.
- A C++ compiler with a makefile would be enough, see “Makefile”.

### Usage

- use Mathematica notebooks to generate codes in ./deps-mma
- copy the generated .txt files into src/aPWD3_part*.cpp
- set N_channels in src/main_part*.cpp
- compile and run
- save channel information and (2J,P,2T) in a txt file

### Method

We follow the "aPWD" method developed by J. Golak et.al.
See "A new way to perform partial wave decompositions of few-nucleon forces" if interested.
We specially thank Professor Kacper Topolnicki for discussions.

### Code Structure

- src/*: C++ codes.
- deps-mma/: Mathematica notebooks to generate spin-projected matrix elements.
- tool/*.py: some simple python tools to generate mesh points, channels index and plot.
- Makefile: you may use it to compile.
- data/: directory of storing 3BMEs.
- script/*.sh: scripts to submit the job on supercomputer.
- docs/Theory.md: details of the aPWD theory.

### Citation

If you benifite from this code during research, please cite as:

- Rongzhe Hu. Three-Body-Jacobi: A code for calculating chiral 3NF matrix element in Jacobi coordinates. [https://github.com/hhurongzhe/Three-Body-Jacobi]

### Remarks

This code is mainly designed for generating non-locally regulated 3N interaction in the Jacobi coordinate (up to N2LO). However, it is possible to extend for local 3N interaction or for N3LO terms.

If you have any needs or questions, just contact me: rongzhe_hu@pku.edu.cn or rongzhehuu@gmail.com.
