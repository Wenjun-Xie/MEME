#!/usr/bin/bash

module load mpich

mpif90 -c    global.f90
mpif90 -c    hamiltonian.f90
mpif90 -c    mc_sampling.f90
mpif90 -c    ensemble_average.f90

mpif90 -o main   *.o main.f90
mpif90 -o energy   *.o energy.f90
rm *.o *.mod

mpirun -np 16 ./main
