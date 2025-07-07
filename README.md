# Three-Body-Jacobi

A high-performance C++17 implementation for calculating the matrix elements of chiral three-nucleon (3N) interactions within an LS-coupled three-body partial-wave basis in Jacobi coordinates.

## Background

The calculation is based on the "Automated-Partial-Wave-Decomposition" (aPWD) method. This approach provides an efficient way to handle the complex spin-isospin algebra and angular momentum coupling involved in few-body force calculations.

For a detailed theoretical description, please refer to:

- **Our documentation**: [`docs/Theory.md`](https://www.google.com/search?q=docs/Theory.md)
- **The original paper**: J. Golak, et al., "A new way to perform partial wave decompositions of few-nucleon forces," _Eur. Phys. J. A_ **43**, _241–250 (2010)_.

## Prerequisites

To build and run this code, you will need:

- **C++17 compiler** (e.g., GCC12 or higher)
- **OpenMP** for parallel computation
- **Wolfram Mathematica** (for generating code snippets from the symbolic notebooks)
- **Python 3** (for using the scripts in the `tool/` directory)

## Usage

Follow these steps:

#### Step 1: Clone the Repository

```bash
git clone https://github.com/hhurongzhe/Three-Body-Jacobi.git
cd Three-Body-Jacobi
```

#### Step 2: Generate C++ Code Snippets (Mathematica)

The core expressions for the partial-wave-projected matrix elements are generated symbolically.

1.  Navigate to the `deps-mma/` directory.
2.  Run the Mathematica notebooks (e.g., `aPWD3_c1.nb`, `aPWD3_c3.nb`, etc.). You need to edit the path in them before run.
3.  This will produce several `.txt` files (e.g., `apwd_Gt_c1_twoJ1_P1_twoT1.txt`), which contain the C++ code for the calculated expressions.

#### Step 3: Integrate Generated Code

The generated code snippets must be manually inserted into the C++ source files.

- For each generated `.txt` file in `deps-mma/`, copy its contents.
- Paste the contents into the corresponding C++ implementation file in `src/`. For example, content from `apwd_Gt_c1_twoJ1_P1_twoT1.txt` should be placed within `src/aPWD3_part_c1.cpp`.

#### Step 4: Configure Calculation Parameters

Before compiling, you need to set the number of channels for your calculation.

- Open the relevant main program file, e.g., `src/main_part_c1.cpp`.
- Set the `N_channels` variable to your desired value.

#### Step 5: Compile the Code

The provided `Makefile` is configured to build the executables.

- In the root directory of the project, run `make`. This will compile the source code and create executables for each part in `build/` directory, `build/apwd3-c1.x`, `build/apwd3-c3.x`, `build/apwd3-c4.x`, `build/apwd3-cD.x`, `build/apwd3-cE.x` and `build/apwd3-benchmark.x`.

```bash
# To build all executables
make -j

# To do a benchmark calculation
./build/apwd3-benchmark.x
```

#### Step 6: Run the Calculation

Execute the compiled binary to start the calculation. The output data will be stored in the `data/` directory. Remember to also save the channel information and the quantum numbers (2J, P, 2T) that you used.

```bash
# Example of calculating c1 part
./build/apwd3-c1.x
```

The `script/` directory contains example shell scripts for submitting jobs to a high-performance computing cluster.

## Code Structure

The project is organized into the following directories:

```
.
├── src/            # C++ source code (.cpp, .hpp)
├── deps-mma/       # Mathematica notebooks for symbolic code generation
├── data/           # Default directory for output data (3BMEs)
├── tool/           # Python helper scripts (mesh generation, plotting, etc.)
├── script/         # Example job submission scripts for HPC systems
├── docs/           # Detailed theoretical documentation
├── Makefile        # Build script for compiling the project
└── README.md       # This file
```

## How to Cite

If you use this code in your research, please cite it.

#### Text Citation:

Rongzhe Hu. _Three-Body-Jacobi: A code for calculating chiral 3NF matrix elements in Jacobi coordinates_. [https://github.com/hhurongzhe/Three-Body-Jacobi](https://github.com/hhurongzhe/Three-Body-Jacobi).

#### BibTeX Entry:

```bibtex
@software{Hu2025ThreeBodyJacobi,
  author       = {Hu, Rongzhe},
  title        = {{Three-Body-Jacobi: A code for calculating chiral 3NF matrix elements in Jacobi coordinates}},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/hhurongzhe/Three-Body-Jacobi}
}
```

## License

This project is licensed under the terms of MIT LICENSE.

## Acknowledgments

We extend special thanks to Professor Kacper Topolnicki for valuable discussions.

## Contact

For questions, bug reports, or further interests, please feel free to contact Rongzhe Hu at:

- **Email**: `rongzhe_hu@pku.edu.cn` or `rongzhehuu@gmail.com`
