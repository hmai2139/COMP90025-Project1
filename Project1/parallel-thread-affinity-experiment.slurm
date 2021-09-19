#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=parallel%j.out
#SBATCH --partition=physical
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=parallel

module load gcc/8.3.0
g++ -fopenmp -Wall -O3 -o Project1 Project1.cpp

echo -e "OMP_PROC_BIND=master\n"

echo - e "1 thread(s)\n"
OMP_PROC_BIND=master OMP_NUM_THREADS = 1 ./Project1 < mseq-big13-example.dat

echo - e "2 thread(s)\n"
OMP_PROC_BIND=master OMP_NUM_THREADS = 2 ./Project1 < mseq-big13-example.dat

echo - e "4 thread(s)\n"
OMP_PROC_BIND=master OMP_NUM_THREADS = 4 ./Project1 < mseq-big13-example.dat

echo - e "8 thread(s)\n"
OMP_PROC_BIND=master OMP_NUM_THREADS = 8 ./Project1 < mseq-big13-example.dat

echo -e "---------------------\n"

echo -e "OMP_PROC_BIND=spread\n"

echo - e "1 thread(s)\n"
OMP_PROC_BIND=spread OMP_NUM_THREADS = 1 ./Project1 < mseq-big13-example.dat

echo - e "2 thread(s)\n"
OMP_PROC_BIND=spread OMP_NUM_THREADS = 2 ./Project1 < mseq-big13-example.dat

echo - e "4 thread(s)\n"
OMP_PROC_BIND=spread OMP_NUM_THREADS = 4 ./Project1 < mseq-big13-example.dat

echo - e "8 thread(s)\n"
OMP_PROC_BIND=spread OMP_NUM_THREADS = 8 ./Project1 < mseq-big13-example.dat

echo -e "---------------------\n"

echo -e "OMP_PROC_BIND=close\n"

echo - e "1 thread(s)\n"
OMP_PROC_BIND=close OMP_NUM_THREADS = 1 ./Project1 < mseq-big13-example.dat

echo - e "2 thread(s)\n"
OMP_PROC_BIND=close OMP_NUM_THREADS = 2 ./Project1 < mseq-big13-example.dat

echo - e "4 thread(s)\n"
OMP_PROC_BIND=close OMP_NUM_THREADS = 4 ./Project1 < mseq-big13-example.dat

echo - e "8 thread(s)\n"
OMP_PROC_BIND=close OMP_NUM_THREADS = 8 ./Project1 < mseq-big13-example.dat

echo -e "---------------------\n"