#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=parallel%j.out
#SBATCH --partition=physical
#SBATCH --time=0:20:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=parallel

module load gcc/8.3.0
g++ -fopenmp -lnuma -Wall -Ofast -o Project1 Project1.cpp
export OMP_PLACES=cores
export OMP_PROC_BIND=close
cat mseq-big13-example.dat | ./Project1
