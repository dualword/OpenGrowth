#!/bin/bash
#SBATCH -n 16
#SBATCH -t 07-00:00:00 
#SBATCH -p general
#SBATCH --mem=20000MB
#SBATCH -e LIGANDNAME.e
#SBATCH -o LIGANDNAME.o
#SBATCH -J LIGANDNAME

Name=LIGANDNAME

module load hpc/gaussian-09
g09 ${Name}.com ${Name}.log

