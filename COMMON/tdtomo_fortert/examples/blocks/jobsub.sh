#!/bin/sh

# Grid Engine Options
#$ -N mksamples
#$ -cwd
#$ -l h_rt=16:00:00
#$ -l h_vmem=2G
#$ -pe mpi 32
#$ -R y

# Initialise the environment
#. /etc/profile
#. /etc/profile.d/modules.sh

# Load the intel and MPI modules
#module load intel
#module load openmpi/1.6.5
#module load igmm/mpi/gcc/mpich/3.1.4

# Run executable
mpiexec -np $NSLOTS ./mksamples 1 
