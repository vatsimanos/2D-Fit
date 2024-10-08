#!/bin/bash -l
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=72
#SBATCH --partition=multinode
#SBATCH --time=02:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
module load openmpi/4.1.2-gcc11.2.0


for ((K=1;K<=1;K++))
do

srun ./2DFit64

done

echo "Stop job :"`date`