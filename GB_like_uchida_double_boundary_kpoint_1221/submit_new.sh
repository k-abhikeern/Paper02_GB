#!/bin/sh

mpirun -np 48 /home/abhikeern/softwares/lammps-5Jun19/src/lmp_mpi -in in.attempt_SLG_NMD_GB > out.dat
