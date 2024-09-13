#!/bin/bash
#$ -N graphene 
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24
#$ -q all.q
/opt/openmpi/bin/mpirun  -np $NSLOTS /share/apps/lammps/lmp_mpi -in  in.attempt_SLG_NMD_GB >  out.dat
