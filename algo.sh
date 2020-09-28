#!/bin/bash

#SBATCH --job-name=ElEnvi
#SBATCH --partition=lr3
#SBATCH --qos=condo_axl
#SBATCH --account=lr_axl
#SBATCH --nodes=1
#SBATCH --time=00:30:00

icc -o probando.exe probando.c -mkl
time mpirun -n 1   ./probando.exe

