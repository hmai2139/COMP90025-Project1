#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=sequential%j.out
#SBATCH --partition=physical
#SBATCH --time=0:40:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=sequential

module load gcc/8.3.0
g++ -fopenmp -Wall -O3 -o kseqalign kseqalign.cpp
./kseqalign < mseq-big13-example.dat
