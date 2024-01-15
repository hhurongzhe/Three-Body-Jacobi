# Three-Body-Jacobi

Calculating matrix elements of chiral 3N interactions under JJ-coupled three-body partial-wave basis ( in Jacobi coordinates ), written in C++.



## Requirements

A C++ compiler with a makefile would be enough.

However, I use Xmake, which is a lightweight and fast building utility based on lua.



## Method

To do partial-wave projection for 3N interactions, we follow the "aPWD" method developed by J. Golak ( see "A new way to perform partial wave decompositions of few-nucleon forces" ). Unlike the description in the original paper, we perform partial-wave decomposition under JJ-coupled 3N states directly ( usually denoted by $\ket{\alpha}$ ). But the formalism is almost the same.



## Structure of Code

- src/constants.hpp: definitions of LECs, mesh points and other things. You may need to change these for your use.
- src/precalculate.hpp: some functions and structures used to precalculate and store in the beginning of calculation. This is necessary for speeding up the code.
- src/aPWD3.hpp: doing partial-wave decomposition for 3N force under JJ-coupled states. This file is rather big because we generate it by Mathematica.
- src/main.hpp: the main function, doing angular integration and calculating in specific channels.
- tools/*.py: some simple python tools to generate mesh points, channels index and some other things useful.
- apwd3.py: a simple and direct code to do above things using python. This is just a template code for check, running very slow.
- xmake.lua: Xmake’s project description file.



## Remarks

This code is mainly designed for generating non-locally regulated 3N interaction in the Jacobi coordinate (up to N2LO). However, it is possible to extend for local 3N interaction or for N3LO terms.

If you have any needs or questions, just contact me: rongzhe_hu@pku.edu.cn


