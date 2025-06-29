#!/bin/bash
#SBATCH --partition c128m512
#SBATCH -o ./output/output-%j-%N.out
#SBATCH -J apwd3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64

export OMP_NUM_THREADS=48
export OMP_STACKSIZE=1g
ulimit -s unlimited


part="cE"

./build/linux/x86_64/release/apwd3-$part.x  >  output/log-xmake-$part.out


