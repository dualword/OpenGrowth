#!/bin/bash

#This script converts a file in .pdb to .xyz, creates an input file
#for Gaussian and the submission script, and submits the latter.
#OpenBabel needs to be installed and available.

module load hpc/openbabel-2.3.1

for Protein in 1HPV
do
##########################
File=${Protein}_Ligand
Charge=` grep "${Protein}" Charges.txt | awk '{ print $3 ; }'`

obabel ${File}.pdb -O ${File}.xyz

#Create the Gaussian input file.
cat << EOF > ${File}.com
%Nprocshared=16
%mem=16000MB
#P M062X/6-31+G** Opt SCF=(MaxCycle=500) SCRF=(IEFPCM,Solvent=Water)
# Test Units(Ang,Deg)

Input

${Charge}       1
EOF

sed '1,2d' ${File}.xyz >> ${File}.com

#Create the script and clean.
sed "s/LIGANDNAME/${File}/g" 00-Files/Run_01-Opt.sh > Run_01-Opt_${File}.sh
chmod u+x Run_01-Opt_${File}.sh
rm ${File}.xyz
sleep 1s

sbatch Run_01-Opt_${File}.sh
##########################
done

