#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=parallel%j.out
#SBATCH --partition=physical
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=parallel

module load gcc/8.3.0
g++ -fopenmp -Wall -O3 -o Project1 Project1.cpp

echo -e "OMP_PROC_BIND=master\n"
OMP_PROC_BIND=master ./Project1 < mseq-big13-example.dat
echo -e "---------------------\n"

echo -e "OMP_PROC_BIND=spread\n"
OMP_PROC_BIND=spread ./Project1 < mseq-big13-example.dat
echo -e "---------------------\n"

echo -e "OMP_PROC_BIND=close\n"
OMP_PROC_BIND=close ./Project1 < mseq-big13-example.dat
echo -e "---------------------\n"