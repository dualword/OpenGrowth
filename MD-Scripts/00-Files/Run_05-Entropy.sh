#!/bin/bash
#SBATCH -n 20
#SBATCH -t 03-00:00:00 
#SBATCH -p general
#SBATCH --mem=50000MB
#SBATCH --array=1-50
#SBATCH -J LIGANDNAME_Entropy
#SBATCH -e LIGANDNAME_%a_Entropy.e
#SBATCH -o LIGANDNAME_%a_Entropy.o

module load centos6/gcc-4.8.2
module load centos6/openmpi-1.6.5_gcc-4.8.0
module load centos6/amber12_openmpi-1.6.5_intel-13.0.079
AMBERHOME=/n/sw/centos6/amber12_openmpi-1.6.5_intel-13.0.079

Protein=PROTEINNAME
Ligand=LIGANDNAME
LigandCharge=CHARGEVALUE
Iteration=${SLURM_ARRAY_TASK_ID}
File=${Ligand}
Cores=20

cd ${File}_${Iteration}/MMPBSA
#Clean everything just in case.
rm -rf Mdcrd_${File}_${Iteration}/ Files_${File}_${Iteration}_Entropy/ MMPBSA_${File}_${Iteration}_Entropy.dat Files_${File}_${Iteration}_Entropy.tar.gz _MMPBSA*

#The coordinates have been compressed previously so uncompress them.
tar -zvxf Mdcrd_${File}_${Iteration}.tar.gz

#Run the entropy calculations.
mpirun -np ${Cores} MMPBSA.py.MPI -O -i ../../00-Files/MMPBSA_Entropy.in -o MMPBSA_${File}_${Iteration}_Entropy.dat -cp Complex.prmtop -rp Receptor.prmtop -lp Ligand.prmtop -y Mdcrd_${File}_${Iteration}/*.mdcrd
sleep 10s

#Clean.
mkdir Files_${File}_${Iteration}_Entropy
mv _MMPBSA* Files_${File}_${Iteration}_Entropy/.
tar czfv Files_${File}_${Iteration}_Entropy.tar.gz Files_${File}_${Iteration}_Entropy/
rm -rf Files_${File}_${Iteration}_Entropy/
rm -rf Mdcrd_${File}_${Iteration}/

