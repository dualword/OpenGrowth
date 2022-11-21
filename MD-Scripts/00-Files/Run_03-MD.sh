#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --gres=gpu
#SBATCH -t 02-00:00:00 
#SBATCH -p gpuqueue
#SBATCH --mem=10000MB
#SBATCH --array=1-50
#SBATCH -J LIGANDNAME
#SBATCH -e LIGANDNAME_%a.e
#SBATCH -o LIGANDNAME_%a.o

module load centos6/gcc-4.8.2
source /n/home04/ncheron/Programs/Gromacs_5.0.2/Threadmpi/bin/GMXRC
#This Gromacs version was installed with the following cmake command:
#cmake .. -DGMX_MPI=OFF -DGMX_GPU=ON -DGMX_THREAD_MPI=ON -DGMX_OPENMP=ON -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/n/home04/ncheron/Programs/Gromacs_5.0.2/Threadmpi/

Protein=PROTEINNAME
Ligand=LIGANDNAME
LigandCharge=CHARGEVALUE
Iteration=${SLURM_ARRAY_TASK_ID}
Extension=
OMP_NUM_THREADS=16
File=${Ligand}
State=STATE

#We perform a loop as long as the simulation didn't succeed.
Success=` grep "Finished mdrun on rank" ${File}_${Iteration}/08-FullMD.log | wc -l `
while [ ${Success} == 0 ]; do
####################################

#Clean everything.
rm -rf ${File}_${Iteration}
echo "*******************" >> ${File}_${Iteration}.e
echo "**** NEW START ****" >> ${File}_${Iteration}.e
echo "*******************" >> ${File}_${Iteration}.e

mkdir ${File}_${Iteration}
cp 00-Files/*.mdp ${File}_${Iteration}/.

#Get the protein structure.
cp 00-Files/Structures/${Protein}_Protein.pdb ${File}_${Iteration}/.

#Define the protonation state. This is specific to this protein (HIV1
#protease). Two aspartic acid residues are in the active site and four
#protonation states can be observed. State=1 means that both are
#unprotonated, State=2 means only ASP25 from chain A is protonated,
#State=3 means only ASP25 from chain B and State=4 means both are protonated.
#if [ ${State}==1 ]; then we use ASP/ASP.
if [ ${State} == 2 ]; then
	sed -i "s/ASP A  25/ASH A  25/g" ${File}_${Iteration}/${Protein}_Protein.pdb
elif [ ${State} == 3 ]; then
	sed -i "s/ASP B  25/ASH B  25/g" ${File}_${Iteration}/${Protein}_Protein.pdb
elif [ ${State} == 4 ]; then
	sed -i "s/ASP A  25/ASH A  25/g" ${File}_${Iteration}/${Protein}_Protein.pdb
	sed -i "s/ASP B  25/ASH B  25/g" ${File}_${Iteration}/${Protein}_Protein.pdb
fi

#Get the water molecules.
cp 00-Files/Structures/${Protein}_Water.pdb ${File}_${Iteration}/.
sed -i -e "s/O   HOH/OW  SOL/g" -e "s/H1  HOH/HW1 SOL/g" -e "s/H2  HOH/HW2 SOL/g" ${File}_${Iteration}/${Protein}_Water.pdb

#Get the ligand files.
cp 00-Files/Structures/${Protein}_Ligand.pdb ${File}_${Iteration}/${Ligand}.pdb
cp 00-Files/Structures/${Ligand}_GMX.itp ${File}_${Iteration}/.
cd ${File}_${Iteration}

echo "#####################"
echo " "
echo "Creation of Topology"
echo " "
echo "#####################"
#We use the Amber03 force field and the TIP3P water model.
gmx${Extension} pdb2gmx -f ${Protein}.pdb -p ${File}.top -o 01-${Protein}.pdb -ignh << EOI
1
1
EOI
sleep 5s

#Cosmetic change if there is only one chain in the protein.
Chain_Number=` ls ${Protein}_Protein_chain_*.itp | wc -l `
if [ ${Chain_Number} == 0 ]; then
	mv posre.itp 01-posre.itp
	sed -i "s/posre.itp/01-posre.itp/g" ${File}.top
fi

echo "###########################################"
echo " "
echo "Prepare and clean the new pdb and top files"
echo " "
echo "###########################################"
sed -i "s/ATOM  /HETATM/g" ${Ligand}.pdb
cat 01-${Protein}.pdb ${Ligand}.pdb ${Protein}_Water.pdb > 01-Complex.pdb
sed -i -e '/MODEL/d' -e '/END/d' -e '/REMARK/d' -e '/COMPND/d' -e '/AUTHOR/d' -e '/CONECT/d' -e '/MASTER/d' 01-Complex.pdb
mv ${File}.top ${File}-Original.top
cat ${File}-Original.top | sed '/forcefield\.itp\"/a\#include "LIGANDNAME_GMX.itp"' >| ${File}.top
sed -i "s/LIGANDNAME/${Ligand}/g" ${File}.top
echo "${Ligand}         1" >> ${File}.top

echo "##################"
echo " "
echo "Solvation in water"
echo " "
echo "##################"
#Solvation in a dodecahedron of 1.1 nm (11 Angstrom)
gmx${Extension} editconf -f 01-Complex.pdb -c -o 02-${File}.gro -d 1.1 -bt dodecahedron
sleep 5s
gmx${Extension} solvate -cp 02-${File}.gro -cs spc216.gro -p ${File}.top -o 02-${File}_Solvated.gro
sleep 5s

echo "######################"
echo " "
echo "Find the system charge"
echo " "
echo "######################"
if [ ${Chain_Number} == 0 ]
then
	grep "qtot" ${File}-Original.top > Temp.txt
	Temp_Size=` wc -l Temp.txt | awk '{ print $1 ; }' `
	Charge_Array=(` awk '{ print $11 ; }' Temp.txt `)
	Index=`expr ${Temp_Size} - 1`
	ChargeTot=` echo "${Charge_Array[${Index}]} ${LigandCharge}" | awk '{print $1+$2}'`
	rm Temp.txt
else
	ChargeTot=${LigandCharge}
	x=0
	ITP_Array=(` ls ${Protein}_Protein_chain_*.itp `)
	while [ $x -lt ${Chain_Number} ]
	do
		grep "qtot" ${ITP_Array[$x]} > Temp-$x.txt
		Temp_Size=` wc -l Temp-$x.txt | awk '{ print $1 ; }' `
		Charge_Array=(` awk '{ print $11 ; }' Temp-$x.txt `)
		Index=`expr ${Temp_Size} - 1`
		ChargeTotNew=` echo "${Charge_Array[${Index}]} ${ChargeTot}" | awk '{print $1+$2}'`
		ChargeTot=${ChargeTotNew}
		rm Temp-$x.txt
		x=`expr $x + 1`
	done
fi

if [ ${ChargeTot} -gt 0 ]
then
        NumberPositive=0
        NumberNegative=${ChargeTot}
elif [ ${ChargeTot} -lt 0 ]
then
        NumberPositive=` echo "scale=1;((-1)*${ChargeTot})" | bc `
        NumberNegative=0
else
        NumberPositive=0
        NumberNegative=0
fi
echo "***** Ligand Charge  = ${LigandCharge} *****"
echo "***** Total Charge   = ${ChargeTot} *****"
echo "***** Number of Positive Ions = ${NumberPositive} *****"
echo "***** Number of Negative Ions = ${NumberNegative} *****"

echo "################"
echo " "
echo "Addition of ions"
echo " "
echo "################"
touch 03-Ions.mdp
gmx${Extension} grompp -v -f 03-Ions.mdp -c 02-${File}_Solvated.gro -p ${File}.top -o 03-Ions.tpr -maxwarn 1	#non-matching atom names
sleep 5s
gmx${Extension} genion -s 03-Ions.tpr -o 03-${File}_Solvated-Ions.gro -p ${File}.top -pname NA -nname CL -np ${NumberPositive} -nn ${NumberNegative} << EOI
SOL
EOI
sleep 5s
mv mdout.mdp 03-Ions_mdout.mdp

if [ ${ChargeTot} == 0 ]; then
	sed -i "s/Water_and_ions/SOL           /g" *.mdp
fi

echo "#########################################"
echo " "
echo "Energy minimization with steepest descent"
echo " "
echo "#########################################"
touch 04-EM-1-Steep.mdp
gmx${Extension} grompp -v -f 04-EM-1-Steep.mdp -c 03-${File}_Solvated-Ions.gro -o 04-EM-1-Steep.tpr -p ${File}.top
sleep 5s
gmx${Extension} mdrun -v -deffnm 04-EM-1-Steep -nb gpu
sleep 5s
mv mdout.mdp 04-EM-1-Steep_mdout.mdp

echo "###########################################"
echo " "
echo "Energy minimization with conjugate gradient"
echo " "
echo "###########################################"
touch 04-EM-2-CG.mdp
gmx${Extension} grompp -v -f 04-EM-2-CG.mdp -c 04-EM-1-Steep.gro -o 04-EM-2-CG.tpr -p ${File}.top
sleep 5s
gmx${Extension} mdrun -v -deffnm 04-EM-2-CG -nb gpu
sleep 5s
mv mdout.mdp 04-EM-2-CG_mdout.mdp

echo "##############"
echo " "
echo "Analysis 04-EM"
echo " "
echo "##############"
gmx${Extension} energy -f 04-EM-1-Steep.edr -o 04-EM-1-Steep_Potential.xvg << EOI
Potential
EOI
gmx${Extension} energy -f 04-EM-2-CG.edr -o 04-EM-2-CG_Potential.xvg << EOI
Potential
EOI
gmx${Extension} editconf -f 04-EM-2-CG.gro -o 04-EM-2-CG.pdb
### xmgrace 04-EM-1-Steep_Potential.xvg
### xmgrace 04-EM-2-CG_Potential.xvg

echo "###################"
echo " "
echo "Generate index file"
echo " "
echo "###################"
gmx${Extension} make_ndx -f 04-EM-2-CG.gro -o 04-index.ndx << EOI
1 | 13
q
EOI

echo "#########################################"
echo " "
echo "Short NVT MD run with position restraints"
echo " "
echo "#########################################"
touch 05-NVT-PR.mdp
gmx${Extension} grompp -v -f 05-NVT-PR.mdp -n 04-index.ndx -o 05-NVT-PR.tpr -c 04-EM-2-CG.gro -r 04-EM-2-CG.gro -p ${File}.top
sleep 5s
gmx${Extension} mdrun -v -deffnm 05-NVT-PR -nb gpu
sleep 5s
mv mdout.mdp 05-NVT-PR_mdout.mdp

echo "##################"
echo " "
echo "Analysis 05-NVT-PR"
echo " "
echo "##################"
gmx${Extension} energy -f 05-NVT-PR.edr -o 05-NVT-PR_Temperature.xvg << EOI
Temperature
EOI
### xmgrace 05-NVT-PR_Temperature.xvg

echo "##########################################"
echo " "
echo "NVT MD with simulated annealing without PR"
echo " "
echo "##########################################"
touch 06-NVT-SA.mdp
gmx${Extension} grompp -v -f 06-NVT-SA.mdp -n 04-index.ndx -o 06-NVT-SA.tpr -c 05-NVT-PR.gro -t 05-NVT-PR.cpt -p ${File}.top
sleep 5s
gmx${Extension} mdrun -v -deffnm 06-NVT-SA -nb gpu
sleep 5s
mv mdout.mdp 06-NVT-SA_mdout.mdp

echo "##################"
echo " "
echo "Analysis 06-NVT-SA"
echo " "
echo "##################"
gmx${Extension} energy -f 06-NVT-SA.edr -o 06-NVT-SA_Temperature.xvg << EOI
Temperature
EOI
### xmgrace 06-NVT-SA_Temperature.xvg

echo "##################################"
echo " "
echo "Short NPT MD run for equilibration"
echo " "
echo "##################################"
touch 07-NPT.mdp
gmx${Extension} grompp -v -f 07-NPT.mdp -n 04-index.ndx -o 07-NPT.tpr -c 06-NVT-SA.gro -t 06-NVT-SA.cpt -p ${File}.top -maxwarn 1	#Using Berendsen pressure coupling invalidates the true ensemble for the thermostat
sleep 5s
gmx${Extension} mdrun -v -deffnm 07-NPT -nb gpu
sleep 5s
mv mdout.mdp 07-NPT_mdout.mdp

echo "###############"
echo " "
echo "Analysis 07-NPT"
echo " "
echo "###############"
gmx${Extension} energy -f 07-NPT.edr -o 07-NPT_Temperature.xvg << EOI
Temperature
EOI
gmx${Extension} energy -f 07-NPT.edr -o 07-NPT_Pressure.xvg << EOI
Pressure
EOI
gmx${Extension} energy -f 07-NPT.edr -o 07-NPT_Density.xvg << EOI
Density
EOI
### xmgrace 07-NPT_Temperature.xvg
### xmgrace 07-NPT_Pressure.xvg
### xmgrace 07-NPT_Density.xvg

echo "###############"
echo " "
echo "Perform Full MD"
echo " "
echo "###############"
touch 08-FullMD.mdp
gmx${Extension} grompp -v -f 08-FullMD.mdp -n 04-index.ndx -o 08-FullMD.tpr -c 07-NPT.gro -t 07-NPT.cpt -p ${File}.top
sleep 5s
gmx${Extension} mdrun -v -deffnm 08-FullMD -nb gpu
sleep 5s
mv mdout.mdp 08-FullMD_mdout.mdp

echo "#####################"
echo " "
echo "Compress output files"
echo " "
echo "#####################"
#If for some reasons the following fails, you can try with "-s 05-NVT-PR.tpr" for example.
gmx${Extension} trjconv -s 02-${File}.gro -n 04-index.ndx -f 05-NVT-PR.trr -o 05-NVT-PR.xtc -pbc nojump -center << EOI
13
Protein_LIG
EOI
gmx${Extension} trjconv -s 02-${File}.gro -n 04-index.ndx -f 06-NVT-SA.trr -o 06-NVT-SA.xtc -pbc nojump -center << EOI
13
Protein_LIG
EOI
gmx${Extension} trjconv -s 02-${File}.gro -n 04-index.ndx -f 07-NPT.trr -o 07-NPT.xtc -pbc nojump -center << EOI
13
Protein_LIG
EOI
#http://www.gromacs.org/Documentation/Terminology/Periodic_Boundary_Conditions
gmx${Extension} trjconv -s 02-${File}.gro -n 04-index.ndx -f 08-FullMD.trr -o 08-FullMD_Center1.xtc -pbc nojump << EOI
System
EOI
gmx${Extension} trjconv -s 02-${File}.gro -n 04-index.ndx -f 08-FullMD_Center1.xtc -o 08-FullMD.xtc -pbc whole -center << EOI
13
Protein_LIG
EOI
rm 05-NVT-PR.trr 06-NVT-SA.trr 07-NPT.trr 08-FullMD.trr 08-FullMD_Center1.xtc

echo "############"
echo " "
echo "Make new tpr"
echo " "
echo "############"
gmx${Extension} tpbconv -s 08-FullMD.tpr -n 04-index.ndx -o 08-FullMD_${Protein}_${Ligand}.tpr << EOI
Protein_LIG
EOI

echo "##################"
echo " "
echo "Analysis 08-FullMD"
echo " "
echo "##################"
### Extract values
gmx${Extension} energy -f 08-FullMD.edr -o 08-FullMD_Potential.xvg << EOI
Potential
EOI
gmx${Extension} energy -f 08-FullMD.edr -o 08-FullMD_TotalEnergy.xvg << EOI
Total-Energy
EOI
gmx${Extension} energy -f 08-FullMD.edr -o 08-FullMD_Temperature.xvg << EOI
Temperature
EOI
gmx${Extension} energy -f 08-FullMD.edr -o 08-FullMD_Pressure.xvg << EOI
Pressure
EOI
gmx${Extension} energy -f 08-FullMD.edr -o 08-FullMD_Density.xvg << EOI
Density
EOI
### xmgrace 08-FullMD_Potential.xvg
### xmgrace 08-FullMD_TotalEnergy.xvg
### xmgrace 08-FullMD_Temperature.xvg
### xmgrace 08-FullMD_Pressure.xvg
### xmgrace 08-FullMD_Density.xvg

### Radius of gyration
gmx${Extension} gyrate -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o 08-FullMD_Gyrate-Protein.xvg << EOI
Protein
EOI
gmx${Extension} gyrate -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o 08-FullMD_Gyrate-MainChain.xvg << EOI
MainChain
EOI

### RMSD from starting structure
gmx${Extension} rms -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o 08-FullMD_RMSD-Protein.xvg << EOI
Protein Protein
EOI
gmx${Extension} rms -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o 08-FullMD_RMSD-MainChain.xvg << EOI
MainChain MainChain
EOI
gmx${Extension} rms -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o 08-FullMD_RMSD-Protein_${Ligand}.xvg << EOI
Protein_LIG Protein_LIG
EOI

### To cluster the protein (least squares fit and RMSD / output)
gmx${Extension} cluster -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -cl ${File}_Clusters.pdb -dist ${File}_Cluster_RMSD-Dist.xvg -o ${File}_Cluster_RMSD-Clust.xpm -b 1000 -e 5000 << EOI
Protein
Protein
EOI
gmx${Extension} xpm2ps -f ${Protein}_Cluster_RMSD-Clust.xpm -o ${Protein}_Cluster_RMSD-Clust.eps -rainbow red
mv cluster.log ${File}_Cluster.log

### Number of interactions between two groups
gmx${Extension} hbond -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -num 08-FullMD_HB-Num_Protein_${Ligand}.xvg -hbm 08-FullMD_HB-Map_Protein_${Ligand}.xpm << EOI
Protein 13
EOI
gmx${Extension} xpm2ps -f 08-FullMD_HB-Map_Protein_${Ligand}.xpm -o 08-FullMD_HB-Map_Protein_${Ligand}.eps -rainbow red

### Get frames at 5ns and at the center of cluster
Value=$(sed '1!d' ${File}_Clusters.pdb)
Time=${Value:19:5}
gmx${Extension} trjconv -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o ${File}_${Iteration}_5000ps.pdb -dump 5000 << EOI
Protein_LIG
EOI
gmx${Extension} trjconv -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -o ${File}_${Iteration}_${Time}ps.pdb -dump ${Time} << EOI
Protein_LIG
EOI

### Check that there is no interactions between periodic images
gmx${Extension} mindist -f 08-FullMD.xtc -s 08-FullMD.tpr -od 08-FullMD_Minimal-Periodic-Distance.xvg -pi << EOI
System
EOI

### Number of interactions between two groups
gmx${Extension} hbond -s 08-FullMD.tpr -n 04-index.ndx -f 08-FullMD.xtc -num 08-FullMD_HB-Num_Protein_${Ligand}.xvg -hbm 08-FullMD_HB-Map_Protein_${Ligand}.xpm << EOI
Protein 13
EOI
gmx${Extension} xpm2ps -f 08-FullMD_HB-Map_Protein_${Ligand}.xpm -o 08-FullMD_HB-Map_Protein_${Ligand}.eps -rainbow red

sleep 30s
cd ..
Success=` grep "Finished mdrun on rank" ${File}_${Iteration}/08-FullMD.log | wc -l `

####################################
done

