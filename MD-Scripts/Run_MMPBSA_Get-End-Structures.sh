#!/bin/bash
#SBATCH -n 1
#SBATCH -t 24:00:00 
#SBATCH -p general
#SBATCH --mem-per-cpu=2500MB
#SBATCH -e Structures.e
#SBATCH -o Structures.o
#SBATCH -J Structures

mkdir Structures
for Protein in 1HPV 
do
	for Iteration in `seq 1 50`
	do
		cp ${Protein}_${Iteration}/${Protein}_${Iteration}_5000ps.pdb Structures/.
	done
done

