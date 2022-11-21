#!/bin/bash

#This script creates the submission script and submits it.

for Protein in 1HPV
do
	Charge=` grep "${Protein}" Charges.txt | awk '{ print $3 ; }'`
	LigandName=` grep "${Protein}" Charges.txt | awk '{ print $1 ; }'`
	sed -e "s/PROTEINNAME/${Protein}/g" -e "s/LIGANDNAME/${LigandName}/g" -e "s/CHARGEVALUE/${Charge}/g" 00-Files/Run_04-MMPBSA.sh > Run_04-MMPBSA_${Protein}.sh
	chmod u+x Run_04-Analysis_${Protein}.sh
	sleep 1s
	sbatch Run_04-Analysis_${Protein}.sh
	sleep 1s
	echo "${Protein} done"
done

