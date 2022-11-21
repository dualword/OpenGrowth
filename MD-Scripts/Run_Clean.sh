#!/bin/bash
#SBATCH -n 1
#SBATCH -t 02-00:00:00 
#SBATCH -p general
#SBATCH --mem=3000MB
#SBATCH -J Clean
#SBATCH -e Clean.e
#SBATCH -o Clean.o

#This scripts cleans the folders.

mkdir Old
for Protein in 1HPV
do
	for Iteration in `seq 1 50`
	do
		if [ -d ${Protein}_${Iteration} ]; then
			LigandName=` grep "${Protein}" Charges.txt | awk '{ print $1 ; }'`
			cd ${Protein}_${Iteration}
			mkdir Old_${Protein}_${Iteration}
			mv 01-* 02-* 03-* 04-* 05-* 06-* 07-* Old_${Protein}_${Iteration}/.
		        mv Old_${Protein}_${Iteration}/04-index.ndx .
		        mv Old_${Protein}_${Iteration}/*.xvg .
			mv MMPBSA/Complex.prmtop MMPBSA/Ligand.prmtop MMPBSA/Receptor.prmtop MMPBSA/*.dat MMPBSA/Files_${Protein}_${Iteration}*.tar.gz .
			mv 08-FullMD.cpt 08-FullMD.edr 08-FullMD.gro 08-FullMD_HB-Map_Protein_${LigandName}.xpm 08-FullMD_mdout.mdp 08-FullMD.mdp 08-FullMD_prev.cpt ${Protein}_Cluster.log ${Protein}_Cluster_RMSD-Clust.xpm *-Original.top *.itp ${Protein}_Protein.pdb ${Protein}_${Iteration}*.e ${Protein}_${Iteration}*.o ${Protein}_Water.pdb Frames_${Protein}_${Iteration}.tar.gz MMPBSA ${LigandName}.pdb Old_${Protein}_${Iteration}/.
			cd ..
			mv ${Protein}_${Iteration}/Old_${Protein}_${Iteration} Old/.
			tar czfv ${Protein}_${Iteration}.tar.gz ${Protein}_${Iteration}/
		fi
	done
done

