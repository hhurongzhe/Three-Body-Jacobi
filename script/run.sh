#!/bin/bash
#SBATCH -A hrz
#SBATCH --partition c128m512
##SBATCH --nodelist node8
#SBATCH -o ./output/output-%j-%N.out
#SBATCH -J apwd3
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

export OMP_NUM_THREADS=48
export OMP_STACKSIZE=1g
ulimit -s unlimited


part="cE"

./build/linux/x86_64/release/apwd3-$part.x  >  output/log-xmake-$part.out


