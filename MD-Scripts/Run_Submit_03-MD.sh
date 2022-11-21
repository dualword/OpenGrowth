#!/bin/bash

#This script creates the submission script and submits it.
#The State value is specific to this protein (the HIV-1 protease),
#see in 00-Files/Run_03-MD.sh for more details.

for Protein in 1HPV
do
	Charge=` grep "${Protein}" Charges.txt | awk '{ print $3 ; }'`
	LigandName=` grep "${Protein}" Charges.txt | awk '{ print $1 ; }'`
	State=3
	sed -e "s/PROTEINNAME/${Protein}/g" -e "s/LIGANDNAME/${LigandName}/g" -e "s/CHARGEVALUE/${Charge}/g" -e "s/STATE/${State}/g" 00-Files/Run_03-MD.sh > Run_03-MD_${Protein}.sh
	chmod u+x Run_03-MD_${Protein}.sh
	sleep 1s
	sbatch Run_03-MD_${Protein}.sh
	sleep 1s
	echo "${Protein} done"
done

