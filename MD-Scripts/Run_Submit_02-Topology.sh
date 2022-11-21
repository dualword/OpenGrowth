#!/bin/bash

#This script needs an optimized geometry of the ligand in the .xyz format.
#It creates the submission script and submits it.

for Protein in 1HPV
do
	File=${Protein}_Ligand
	mkdir ${File}
	mv ${File}.xyz ${File}/.
	Charge=` grep "${Protein}" Charges.txt | awk '{ print $3 ; }'`
	sed -e "s/LIGANDNAME/${File}/g" -e "s/CHARGEVALUE/${Charge}/g" 00-Files/Run_02-Topology.sh > ${File}/Run_02-Topology_${File}.sh
	cd ${File}
	chmod u+x Run_02-Topology_${File}.sh
	sleep 1s
	sbatch Run_02-Topology_${File}.sh
	cd ..
done

