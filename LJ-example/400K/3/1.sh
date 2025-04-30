#!/bin/bash

# NOTE: Lines starting with "#SBATCH" are valid SLURM commands or statements,
#       while those starting with "#" and "##SBATCH" are comments.  Uncomment
#       "##SBATCH" line means to remove one # and start with #SBATCH to be a
#       SLURM command or statement.

#SBATCH -J ionic_crystal	 #Slurm job name

# Set the maximum runtime, uncomment if you need it
#SBATCH -t 24:00:00 #Maximum runtime of 24 hours

# Choose partition (queue) to use. Note: replace <partition_to_use> with the name of partition
#SBATCH -p cpu-share

#SBATCH -N 1 -n 40

# Setup runtime environment if necessary
# For example, setup intel MPI environment

# Go to the job submission directory and run your application

module unload ohpc
module load gnu8/8.3.0
module load openmpi3/3.1.4

mpirun -np 40 /home/wliucj/lammps-2Aug2023/src/lmp_mpi -in EMD_enthalpy.lmp

