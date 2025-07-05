#!/bin/bash
#SBATCH --partition c128m512
#SBATCH -o ./output/output-%j-%N.out
#SBATCH -J apwd3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64

ulimit -s unlimited
export OMP_NUM_THREADS=128
export OMP_STACKSIZE=1g

part="c1"
./build/apwd3-$part.x  >  output/log-apwd3-$part.txt
