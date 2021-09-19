#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=parallel%j.out
#SBATCH --partition=physical
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=parallel

module load gcc/8.3.0
g++ -fopenmp -lnuma -Wall -O3 -o Project1 Project1.cpp
echo -e "32 threads\n"
for ((n = 0; n < 10; n++)); do OMP_PLACES=cores OMP_PROC_BIND=close OMP_NUM_THREADS=32 cat mseq-big13-example.dat | ./Project1; done 


