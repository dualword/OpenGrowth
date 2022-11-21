#include "OpenGrowth.h"

void   StoreNames(Parameters const & parameters, string snapshotNames[], string const typeOfSnapshots);
double Energy(Parameters const & parameters, Protein const & protein, Molecule & molecule);
int    StericClash(Parameters const & parameters, Protein const & protein, Molecule const & molecule);
void   OptimizationPosition(Parameters & parameters, Protein const & protein, Molecule & molecule);
float  OptimizationGeom(Parameters const & parameters, Protein const & protein, Molecule & molecule, int & errorFF, int const & optParameter);

// This function prepares the first fragment when SEED or ENERGY mode is selected. The LIGAND parameter must have been defined.
void  FirstFragmentIfKnown(Parameters & parameters, Protein const protein[], Molecule molecule[], double & averageEnergy, string const typeOfSnapshots, int & errorFF)
{
    // Depending on which kind of snapshots we are using, there are not the same number of them.
    int snapshotNumber=0;
    if      (typeOfSnapshots == "rotamers")   { snapshotNumber = parameters.rotamersNumber;   }
    else if (typeOfSnapshots == "conformers") { snapshotNumber = parameters.conformersNumber; }

    // Store names.
    string ligandNames[snapshotNumber];                                                                            // Array for the ligand names; we assume the same number of ligands than snapshots.
    if      (typeOfSnapshots == "rotamers")   { StoreNames(parameters, ligandNames, "ligandsRotamer");   }         // Store the ligand names for rotamers.
    else if (typeOfSnapshots == "conformers") { StoreNames(parameters, ligandNames, "ligandsConformer"); }         // Store the ligand names for conformers.

    // Array of ligands.
    Molecule ligands[snapshotNumber];

    for (int s=0 ; s < snapshotNumber ; s++) {
        // Create an object for OpenBabel.
        OpenBabel::OBConversion obConversion;
        OpenBabel::OBFormat *format = obConversion.FormatFromExt(ligandNames[s]);
        obConversion.SetInFormat(format);
        obConversion.ReadFile(&ligands[s].obMol, ligandNames[s]);

        // Open the ligand file, and check if it has been opened.
        ifstream ligandFile(ligandNames[s].c_str(), ios::in);
        if(!ligandFile ) {
            cerr << "ERROR: The file " << ligandNames[s] << " is missing, or I can't read it." << endl;
            exit (-1);
            }

        // Prepare the fragment.
        ligands[s].growth.clear();
        ligands[s].fragmentIndex.clear();
        ligands[s].fragmentNeighbours.clear();
        int IsARing=0;
        int IsAromatic=0;

        // Read the fragment.
        unsigned int sizeFragment=0;
        string line;
        int lineCounter=1;
        while(getline(ligandFile, line)) {
            istringstream input(line);
            if (lineCounter==1) {
                string sizeValue="";
                input >> sizeValue;
                sizeFragment=atoi(sizeValue.c_str());
                }
            else if (lineCounter==2) {
                string nameValue="";
                input >> nameValue;
                }
            else {
                string atomTypeValue="";
                string xValue="";
                string yValue="";
                string zValue="";
                string growthValue="";
                input >> atomTypeValue >> xValue >> yValue >> zValue >> growthValue;
                int atomGrowth = atoi(growthValue.c_str());
                // Update the dynamic arrays.
                ligands[s].growth.push_back(atomGrowth);
                ligands[s].fragmentIndex.push_back(1);
                ligands[s].fragmentNeighbours.push_back(0);
                }
            lineCounter++;
            }

        // Check if each atom is aromatic or is in a ring
        for (unsigned int j=1; j<=ligands[s].obMol.NumAtoms() ; j++ ) {
            // Get info from the atom.
            OpenBabel::OBAtom *atom;
            atom = ligands[s].obMol.GetAtom(j);
            if (ligands[s].growth[j-1]==1) {
                IsARing += atom->IsInRing();
                IsAromatic += atom->IsAromatic();
                }
            }

        // Define some additional parameters.
        ligands[s].growingSites.clear();
        for (unsigned int j=1; j<=ligands[s].growth.size() ; j++ ) {
            if(ligands[s].growth[j-1]>=1) { ligands[s].growingSites.push_back(j); }
            }
        ligands[s].numberResidues = 1;
        ligands[s].stericClashes = 0;
        if (IsARing == 0)    { ligands[s].ring = 0; }
        else                 { ligands[s].ring = 1; }
        if (IsAromatic == 0) { ligands[s].aromatic = 0; }
        else                 { ligands[s].aromatic = 1; }

        // Check that all molecule size are consistent, and close the file.
        if ( (ligands[s].obMol.NumAtoms()!=sizeFragment) || (ligands[s].obMol.NumAtoms()!=ligands[s].growth.size()) || (ligands[s].growth.size()!=sizeFragment) ) {
            cerr << "ERROR: some of the arrays describing the ligand have not the same size." << endl;
            exit (-1);
            }
        ligandFile.close();

        // The ligand becomes the molecule, and we optimize the molecule.
        molecule[s] = ligands[s];
        molecule[s].interactionEnergy = Energy(parameters, protein[s], molecule[s]);
        molecule[s].stericClashes = StericClash(parameters, protein[s], molecule[s]);
        if ((parameters.optimizationMode/10)>=2) { OptimizationGeom(parameters, protein[s], molecule[s], errorFF, 1); }
        if ((parameters.optimizationMode%10)>=2) { OptimizationPosition(parameters, protein[s], molecule[s]);         }
        }

    // Compute the average energy of the molecules.
    averageEnergy=0;
    long double sumEnergy=0;
    if      (snapshotNumber==1)                   { averageEnergy=molecule[0].interactionEnergy; }
    else if (parameters.averageType=="BOLTZMANN") {
        // Compute the Boltzmann average energy.
        long double partitionFunction=0;
        for (int s=0 ; s < snapshotNumber ; s++) {
            partitionFunction += exp(-BETA*molecule[s].interactionEnergy);
            sumEnergy += molecule[s].interactionEnergy*exp(-BETA*molecule[s].interactionEnergy);
            }
        averageEnergy = sumEnergy/partitionFunction;
        }
    else if (parameters.averageType=="ARITHMETIC") {
        // Compute the arithmetic average energy.
        for (int s=0 ; s < snapshotNumber ; s++) {
            averageEnergy += molecule[s].interactionEnergy;
            }
        averageEnergy = averageEnergy/snapshotNumber;
        }
    else if (parameters.averageType=="LOWESTSCORE") {
        // Keep only the lowest energy.
        averageEnergy = molecule[0].interactionEnergy; 
        for (int s=1 ; s < snapshotNumber ; s++) {
            if(molecule[s].interactionEnergy<averageEnergy) { averageEnergy = molecule[s].interactionEnergy; }
            }
        }

    // Check for steric clashes in all the molecules.
    int stericClashes=0, clashIntra=0, clashInter=0;
    for (int s=0 ; s < snapshotNumber ; s++) {
        if      ( molecule[s].stericClashes==1 ) { clashInter++;               }
        else if ( molecule[s].stericClashes==2 ) { clashIntra++;               }
        else if ( molecule[s].stericClashes==3 ) { clashInter++; clashIntra++; }
        }
    if      (clashIntra==0 && clashInter==0) { stericClashes=0; }
    else if (clashIntra==0 && clashInter>0)  { stericClashes=1; }
    else if (clashIntra>0  && clashInter==0) { stericClashes=2; }
    else if (clashIntra>0  && clashInter>0)  { stericClashes=3; }

    // Display results.
    cout << endl << "Energy of interaction between the protein and the ligand: " << averageEnergy << endl;

    // If steric clashes are found, display a message.
    if      (stericClashes && (parameters.mode=="ENERGY"))  {
        cout << "Steric clashes have been found between the protein and the ligand." << endl;
        cout << "This should not happen. Check your ligand coordinates or change VDW_SCALE." << endl;
        }
    else if (stericClashes && (parameters.mode=="SEED")) {
        cerr << "ERROR: Steric clashes have been found between the protein and the ligand." << endl;
        cout << "This should not happen. Check your ligand coordinates or change VDW_SCALE." << endl;
        exit (-1); 
        }
}

