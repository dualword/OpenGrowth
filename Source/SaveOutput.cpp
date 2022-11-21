#include "OpenGrowth.h"
#include <iomanip>               // For setprecision
#include <openbabel/parsmart.h>

void  OptimizationPosition(Parameters & parameters, Protein const & protein, Molecule & molecule);
float OptimizationGeom(Parameters const & parameters, Protein const & protein, Molecule & molecule, int & errorFF, int const & optParameter);
float LigandConstraints(Parameters const & parameters, Molecule molecule);

// This function saves a molecule to an .xyz output file. The name is the OUTPUT parameter, then "x_s.xyz" where x is a counter and s the snapshot number.
void SaveOutput(Parameters & parameters, Protein const protein[], Molecule molecule[], int const ligandDescription[], double & averageEnergy, int const & counterOutput, string const threeMerList[], string const typeOfSnapshots, int & errorFF)
{
    // Depending on which kind of snapshots we are using, there are not the same number of them.
    int snapshotNumber=0;
    if      (typeOfSnapshots == "rotamers")   { snapshotNumber = parameters.rotamersNumber;   }
    else if (typeOfSnapshots == "conformers") { snapshotNumber = parameters.conformersNumber; }

    // Open the summary file.
    string outputSummary = parameters.outputName + "_Summary.txt";
    ofstream outputSum(outputSummary.c_str(), ios::out|ios::app);
    if(!outputSum) {
        cerr << "ERROR: I cannot write in the output file " << outputSummary << "." << endl;
        exit (-1);
        }
    // Open the summary SMILES file.
    string outputSummarySMILES = parameters.outputName + "_Summary.smi";
    ofstream outputSumSMILES(outputSummarySMILES.c_str(), ios::out|ios::app);
    if(!outputSumSMILES) {
        cerr << "ERROR: I cannot write in the output file " << outputSummarySMILES << "." << endl;
        exit (-1);
        }

    // Counters for the output file.
    ostringstream counterOutputFlow;
    counterOutputFlow << counterOutput;
    string outputNumber = counterOutputFlow.str();

    // Optimize the ligand.
    for (int s=0 ; s < snapshotNumber; s++) {
        if ((parameters.optimizationMode/10)>=1) { OptimizationGeom(parameters, protein[s], molecule[s], errorFF, 1); }
        if ((parameters.optimizationMode%10)>=1) { OptimizationPosition(parameters, protein[s], molecule[s]);         }
        }

    // Recompute the last energy, since it may have changed after the previous optimizations.
    long double newEnergy=0;
    int snapshotBestValue=0;
    if      (snapshotNumber==1)                   { newEnergy=molecule[0].interactionEnergy; }
    else if (parameters.averageType=="BOLTZMANN") {
        // Compute the Boltzmann average energy of the new molecules.
        long double partitionFunction=0;
        for (int s=0 ; s < snapshotNumber ; s++) {
            partitionFunction += exp(-BETA*molecule[s].interactionEnergy);
            newEnergy += molecule[s].interactionEnergy*exp(-BETA*molecule[s].interactionEnergy);
            }
        newEnergy = newEnergy/partitionFunction;
        }
    else if (parameters.averageType=="ARITHMETIC") {
        // Compute the arithmetic average energy of the new molecules.
        for (int s=0 ; s < snapshotNumber ; s++) {
            newEnergy += molecule[s].interactionEnergy;
            }
        newEnergy = newEnergy/snapshotNumber;
        }
    else if (parameters.averageType=="LOWESTSCORE") {
        // Keep only the lowest energy.
        newEnergy = molecule[0].interactionEnergy;
        for (int s=1 ; s < snapshotNumber ; s++) {
            if(molecule[s].interactionEnergy<newEnergy) { newEnergy = molecule[s].interactionEnergy; snapshotBestValue=s; }
            }
        }
    averageEnergy = newEnergy;

    // Write the ligand for all snapshots, and do it only if parameters.smilesOnly is 0.
    for (int s=0 ; s < snapshotNumber && parameters.smilesOnly==0 ; s++) {
        // Counters for the snapshot index.
        ostringstream counterSnapshot;
        counterSnapshot << s;
        string outputSnapshot = counterSnapshot.str();

        // Full output name.
        string outputFileName = parameters.outputName + "_" + outputNumber + "_" + outputSnapshot + ".xyz";

        OpenBabel::OBConversion obConversion;
        obConversion.SetOutFormat("xyz");
        molecule[s].obMol.SetTitle(outputFileName);
        obConversion.WriteFile(&molecule[s].obMol, outputFileName);
        }

    // Get the SMILES string.
    OpenBabel::OBConversion obConversion;
    obConversion.SetOutFormat("smi");
    obConversion.AddOption("n",OpenBabel::OBConversion::OUTOPTIONS); 
    string smilesString = obConversion.WriteString(&molecule[0].obMol);
    smilesString[smilesString.length()-1]=' ';
    // Write in the SMILES file.
    outputSumSMILES << smilesString << "\t" << parameters.outputName + "_" + outputNumber << "\t" << averageEnergy << "\n" ;

    // Define the description size.
    int descriptionSize=0;
    if      ( parameters.growthMode=="RANDOM" || parameters.growthMode=="BIASED" || parameters.growthMode=="FOG" ) { descriptionSize = 3*molecule[0].numberResidues-2;   }
    else if ( parameters.growthMode=="REGROW")                                                                     { descriptionSize = 3*parameters.maxFragment-2;       }

    // 3Mer-Screen: look for the SMARTS patterns.
    int Counting3Mer=0;
    if (parameters.threeMerFile != "") {
        for (int z=0 ; z < parameters.threeMerSize && !Counting3Mer ; z++) {
            OpenBabel::OBSmartsPattern obSmarts;
            obSmarts.Init(threeMerList[z]);
            if (obSmarts.HasMatch(molecule[0].obMol)) { Counting3Mer++; }
            //obSmarts.Match(molecule[0].obMol);
            //vector<vector<int> > obMaplist;
            //obMaplist = obSmarts.GetUMapList();
            //Counting3Mer += obMaplist.size();
            }
        }

    // Ligand constraints
    float constraints = LigandConstraints(parameters, molecule[0]);
    for (int s=0 ; s < snapshotNumber ; s++) {
        molecule[s].internalEnergy = constraints;
        }

    // Write in the summary file (the first line is for the setprecision below). It is not very beautiful, but I am doing it this way because I am trying to have a
    // formatted output file. I am not using \t because depending on the text editor it is not of the same size.
    outputSum << fixed;
    // Name
    outputSum << parameters.outputName + "_" + outputNumber ;
    if      ( counterOutput < 10 )       { outputSum << "            " ; }  // in OpenGrowth.cpp, we said we want 13 " ". Here we have only 12 because 1 is taken by the index.
    else if ( counterOutput < 100 )      { outputSum << "           " ;  }
    else if ( counterOutput < 1000 )     { outputSum << "          " ;   }
    else if ( counterOutput < 10000 )    { outputSum << "         " ;    }
    else if ( counterOutput < 100000 )   { outputSum << "        " ;     }
    else if ( counterOutput < 1000000 )  { outputSum << "       " ;      }
    else if ( counterOutput < 10000000 ) { outputSum << "      " ;       }
    else                                 { outputSum << "     " ;        }
    outputSum << setprecision(2) ;
    // Score
    outputSum << averageEnergy ;
    if      ( averageEnergy > -10 )  { outputSum << "      " ; }
    else if ( averageEnergy > -100 ) { outputSum << "     " ;  }
    else                             { outputSum << "    " ;   }
    // Snapshot with lowest score
    if (parameters.averageType=="LOWESTSCORE") {
        if      ( snapshotBestValue > 10 ) { outputSum << snapshotBestValue << "           "; }
        else                               { outputSum << snapshotBestValue << "            "; }
        }
    // #Residues
    outputSum << molecule[0].numberResidues ;
    if      ( molecule[0].numberResidues <= 9 )  { outputSum << "             " ; }
    else if ( molecule[0].numberResidues <= 99 ) { outputSum << "            " ;  }
    else                                         { outputSum << "           " ;   }
    // #HeavyAtoms
    outputSum << molecule[0].obMol.NumHvyAtoms() ;
    if      ( molecule[0].obMol.NumHvyAtoms() <= 9 )  { outputSum << "               " ; }
    else if ( molecule[0].obMol.NumHvyAtoms() <= 99 ) { outputSum << "              " ;  }
    else                                              { outputSum << "             " ;   }
    // Molecular weight
    outputSum << molecule[0].obMol.GetMolWt() ;
    if      ( molecule[0].obMol.GetMolWt() < 100 )  { outputSum << "                " ; }
    else if ( molecule[0].obMol.GetMolWt() < 1000 ) { outputSum << "               " ;  }
    else                                            { outputSum << "              " ;   }
    // 3Mer
    if      ((parameters.threeMerFile != "") && !Counting3Mer) { outputSum << "PASS     " ; }
    else if ((parameters.threeMerFile != "") && Counting3Mer)  { outputSum << "FAIL     " ; }
    // Constraints
    outputSum << constraints ;
    if      ( constraints < -100 ) { outputSum << "             " ; }
    else if ( constraints < -10  ) { outputSum << "              " ; }
    else if ( constraints < 0    ) { outputSum << "               " ; }
    else if ( constraints < 10   ) { outputSum << "                " ; }
    else if ( constraints < 100 ) { outputSum << "               " ;  }
    else                          { outputSum << "              " ;  }
    // Description
    if (parameters.writeDescription == 1) {
        int spaceToRemove=0;
        for (int i=0 ; i < descriptionSize ; i++) {
            if (ligandDescription[i]>=10)  { spaceToRemove++; }
            if (ligandDescription[i]>=100) { spaceToRemove++; }
            }
        for (int i=0 ; i < descriptionSize ; i++)                                                                 { outputSum << ligandDescription[i] << " "; }
        for (int i=0 ; i < 3*(3*parameters.maxFragment-2)-2*(3*molecule[0].numberResidues-2)-spaceToRemove ; i++) { outputSum << " ";                         }
        }
    // SMILES
    outputSum << smilesString << "\n" ;

    // Close the summary file and the input file.
    outputSum.close();
    outputSumSMILES.close();
}

