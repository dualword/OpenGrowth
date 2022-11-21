#!/bin/bash
#SBATCH -n 16
#SBATCH -t 04:00:00 
#SBATCH -p general
#SBATCH --mem=16000MB
#SBATCH -e Topology_LIGANDNAME.e
#SBATCH -o Topology_LIGANDNAME.o
#SBATCH -J Topology_LIGANDNAME

#This scripts need an optimized structure in the .xyz format
#It first loads Gaussian to calculate the electrostatic potential,
#then prepares the topology files needed for Gromacs and Amber.

Name=LIGANDNAME
Charge=CHARGEVALUE

#Make Gaussian calculation
cat << EOF > ${Name}.com
%Nprocshared=16
%mem=12000MB
#P HF/6-31G* Pop=MK IOp(6/33=2) IOp(6/41=15) IOp(6/42=15) IOp(6/50=1)
# Test Units(Ang,Deg)

Input

${Charge}       1
EOF

sed '1,2d' ${Name}.xyz >> ${Name}.com

cat << EOF >> ${Name}.com

g09.gesp

g09.gesp

EOF

module load hpc/gaussian-09
g09 ${Name}.com ${Name}.log
echo "Gaussian Done"

#Use antechamber to get the charges in the mol2 file.
module load bio/AmberTools12
antechamber -i ${Name}.log -fi gout -o ${Name}.mol2 -fo mol2 -c resp -eq 2 -nc ${Charge} -s 2 -rn LIG
echo "Antechamber done"

#Clean a little bit
mkdir 01-Gaussian
mv ${Name}.xyz ${Name}.com ${Name}.log 01-Gaussian/.
mkdir 02-Antechamber
mv ANTECHAMBER_AC.AC  ANTECHAMBER_AC.AC0  ANTECHAMBER_BOND_TYPE.AC  ANTECHAMBER_BOND_TYPE.AC0  ANTECHAMBER.ESP  ANTECHAMBER_RESP1.IN  ANTECHAMBER_RESP1.OUT  ANTECHAMBER_RESP2.IN  ANTECHAMBER_RESP2.OUT  ANTECHAMBER_RESP.AC ATOMTYPE.INF qout QOUT esout punch 02-Antechamber/.
cp ${Name}.mol2 02-Antechamber/.

#Use acpype.py to prepare the needed files
module load hpc/python-2.7.3
acpype.py -i ${Name}.mol2 -c user
echo "Acpype done"

#Clean again
mv ${Name}.acpype 03-${Name}.acpype
rm ${Name}.mol2
cp 03-${Name}.acpype/${Name}_GMX.itp .
cp 03-${Name}.acpype/${Name}_AC.frcmod .
cp 03-${Name}.acpype/${Name}.mol2 .
#cp ${Name}_GMX.itp ${Name}_AC.frcmod ${Name}.mol2 ../00-Files/Structures/.


