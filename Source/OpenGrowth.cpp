#include "OpenGrowth.h"

void  Parse(int nbarg, char * argv[], Parameters & parameters);
void  PrepareFragments(Parameters & parameters, Molecule fragment[]);
void  StoreNames(Parameters const & parameters, string snapshotNames[], string const typeOfSnapshots);
void  PrepareProtein(Parameters const & parameters, Protein protein[], string const snapshotNames[], string const typeOfSnapshots);
void  Prepare3Mer(Parameters const & parameters, string threeMerList[]);
void  PrepareProbabilityFiles(Parameters & parameters, double probaFirstFrag[], double probaTransition[][MAX_FRAGMENTS]);
void  PrepareRegrowFile(Parameters const & parameters, int ligandDescription[]);
void  FirstFragment(Parameters & parameters, Protein const protein[], Molecule molecule[], double & averageEnergy, string const typeOfSnapshots, Molecule const fragment[], double const probaFirstFrag[], double randomAxe[], double randomPoint[], int ligandDescription[], int & errorFF);
void  FirstFragmentIfKnown(Parameters & parameters, Protein const protein[], Molecule molecule[], double & averageEnergy, string const typeOfSnapshots, int & errorFF);
void  AddFragment(Parameters & parameters, Protein const protein[], Molecule molecule[], double & averageEnergy, string const typeOfSnapshots, Molecule const fragment[], double const probaFirstFrag[], double const probaTransition[][MAX_FRAGMENTS], int ligandDescription[], int & errorFF); 
void  SaveOutput(Parameters & parameters, Protein const protein[], Molecule molecule[], int const ligandDescription[], double & averageEnergy, int const & counterOutput, string const threeMerList[], string const typeOfSnapshots, int & errorFF);
void  CheckGrowth(Molecule const molecule);

//////////////////////////////////////////////////////////////////////////////////////////

// We define a function which compute the distance between two atoms.
double dist(double const atom1[3], double const atom2[3])
{
    double Distance=sqrt( pow(atom1[0]-atom2[0],2) + pow(atom1[1]-atom2[1],2) + pow(atom1[2]-atom2[2],2) );
    return Distance;
}
// We also define a function which compute the square of the distance to avoid computing the square root when we can.
double dist2(double const atom1[3], double const atom2[3])
{
    double Distance2=pow(atom1[0]-atom2[0],2) + pow(atom1[1]-atom2[1],2) + pow(atom1[2]-atom2[2],2);
    return Distance2;
}

// This function calculates the number of non-empty lines in a file.
int SizeFile(string const fileName)
{
    // Open the file and check if it can be opened.
    // We don't verify if the file has been opened: if it hasn't then its size is zero.
    ifstream listFile(fileName.c_str(), ios::in);

    // If a read line is not empty, increase the counter.
    int fileSize=0;
    string line;
    while(getline(listFile, line)) {
        if( line!="" ) { fileSize++; }
        }
    listFile.close();

    return fileSize;
}

//////////////////////////////////////////////////////////////////////////////////////////

int main(int nbarg, char * argv[])
{
    cout << "************************************************************************" << endl;
    cout << "***************                 OpenGrowth                **************" << endl;
    cout << "***************                    v1.0                   **************" << endl;
    cout << "************************************************************************" << endl;
    cout << "***************           Author: Nicolas Chéron          **************" << endl;
    cout << "***************          Eugene Shakhnovich Group         **************" << endl;
    cout << "***************           Harvard CCB, 2015-2017          **************" << endl;
    cout << "***************            License: GNU GPL v3            **************" << endl;
    cout << "************************************************************************" << endl;
    cout << "***************               Please cite:                **************" << endl;
    cout << "* J. Med. Chem. 2016, pp 4171-4188, DOI: 10.1021/acs.jmedchem.5b00886  *" << endl;
    cout << "* J. Chem. Inf. Model. 2017, pp 584-593, DOI: 10.1021/acs.jcim.6b00610 *" << endl;
    cout << "************************************************************************" << endl;
    cout << endl;

    // Parse the input file.
    Parameters parameters;
    Parse(nbarg, argv, parameters);

    // Initialize random numbers.
    parameters.seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 randomTemp(parameters.seed);
    parameters.random = randomTemp;

    // Get the full path.
    char *full_path = realpath(argv[1], NULL);
    string fullpath = full_path;
    unsigned int positionToErase = fullpath.find_last_of("/");
    fullpath.erase(positionToErase);

    // Create the output directory in case it doesn't exist.
    string outputDir = fullpath + "/" + parameters.outputName;
    positionToErase = outputDir.find_last_of("/");
    outputDir.erase(positionToErase);
    mkdir(outputDir.c_str(), 0777);

    // Prepare fragments and proteins.
    Molecule fragment[MAX_FRAGMENTS];                               // Array of fragments
    Protein proteinRotamers[parameters.rotamersNumber];             // Array of protein rotamers
    Protein proteinConformers[parameters.conformersNumber];         // Array of protein conformers
    PrepareFragments(parameters, fragment);                         // Prepare the fragments (store the address, name, ring information, growth mode of atoms)
    if (parameters.rotamersNumber!=0)   { 
        string rotamersNames[parameters.rotamersNumber];                         // Array for the rotamer names
        StoreNames(parameters, rotamersNames, "rotamers");                       // Store the rotamer names
        PrepareProtein(parameters, proteinRotamers, rotamersNames, "rotamers");  // Prepare the protein rotamers (store coordinates and atom types)
        }
    if (parameters.conformersNumber!=0) { 
        string conformersNames[parameters.conformersNumber];                           // Array for the conformer names
        StoreNames(parameters, conformersNames, "conformers");                         // Store the conformer names
        PrepareProtein(parameters, proteinConformers, conformersNames, "conformers");  // Prepare the protein conformers (store coordinates and atom types)
        }

    // Prepare the forbidden 3-mer list.
    if (parameters.threeMerFile!="") {
        cout << "************************** Prepare 3mer file *************************" << endl;
        parameters.threeMerSize=SizeFile(parameters.threeMerFile);
        }
    else                             { parameters.threeMerSize=0; }
    string threeMerList[parameters.threeMerSize];
    Prepare3Mer(parameters, threeMerList);

    // If the growth mode is FOG, store information from the probability libraries.
    double probaFirstFrag[parameters.fragmentListSize];
    double probaTransition[parameters.fragmentListSize][MAX_FRAGMENTS];
    if (parameters.growthMode=="FOG" || parameters.growthMode=="BIASED") { PrepareProbabilityFiles(parameters, probaFirstFrag, probaTransition); }

    // One of the growth mode is "REGROW". It allows to reconstruct ligands, by specifing the numbers which would have been randomly chosen otherwise.
    // We call it the "ligand description", and we define an array for it. The first number is the index of the first fragment. Then, for every new fragment,
    // 3 numbers are given: the index of the hydrogen in the current ligand, the index of the new fragment, the index of the hydrogen in the chosen fragment.
    // Thus, the size of the array is 3*(number of fragments)-2. Below, if REGROW is chosen we store the values, otherwise we fill the array with "999".
    int  ligandDescription[3*parameters.maxFragment-2];
    if   (parameters.growthMode=="REGROW") { PrepareRegrowFile(parameters, ligandDescription);                                  }
    else                                   { for (int i=0 ; i < 3*parameters.maxFragment-2 ; i++) { ligandDescription[i]=999; } }

    // Prepare a summary file. I am doing it this way (not very beautiful) because I am trying to have a formatted output file. I am not using \t because
    // depending on the text editor it is not of the same size.
    string outputSummary = parameters.outputName + "_Summary.txt";
    ofstream outputSum(outputSummary.c_str(), ios::out|ios::app);
    if(!outputSum) {
        cerr << "ERROR: I cannot write in the output file " << outputSummary << "." << endl;
        exit (-1);
        }

    // Get the size of the output file. If it is 0, write the header. This could be a separate function.
    int sizeOutputSum=SizeFile(outputSummary);
    if (sizeOutputSum==0) {
        outputSum << "#LigandName" ;
        for (int unsigned i=0 ; i < parameters.outputName.size()+1+13-11 ; i++) { outputSum << " "; }   // The +1 is for the "_", 13 is the space we want after, 11 is the length of "#LigandName".
        outputSum << "Score      " ; 
        if (parameters.averageType=="LOWESTSCORE") { outputSum << "Snapshot     "; }
        outputSum << "Fragments     Heavy atoms     Molecular weight     " ;
        if (parameters.threeMerFile != "") { outputSum << "3Mer     " ; }
        outputSum << "LigandConstraints   " ; 
        if (parameters.writeDescription == 1) {
            outputSum << "Description" ;
            for (int i=0 ; i < 3*(3*parameters.maxFragment-2)-11 ; i++) { outputSum << " "; }           // 11 is the length of "Description"
            }
        outputSum << "SMILES" << endl;
        }
    outputSum.close();

    cout << endl << "######################################################################" << endl;
    if (parameters.forbiddenFragments.size()!=0) {
        cout << "The following fragments will not be used: " ;
        for (unsigned int j=0; j<parameters.forbiddenFragments.size() ; j++ ) {
            cout << parameters.forbiddenFragments[j] << " " ;
            }
        cout << endl << "######################################################################" << endl;
        }

    // If we only want the energy.
    if (parameters.mode=="ENERGY") {
        int errorFF=0;
        if (parameters.rotamersNumber!=0) {
            double averageEnergy=0;
            Molecule molecule[parameters.rotamersNumber];
            FirstFragmentIfKnown(parameters, proteinRotamers, molecule, averageEnergy, "rotamers", errorFF);
            }
        if (parameters.conformersNumber!=0) {
            double averageEnergy=0;
            Molecule molecule[parameters.conformersNumber];
            FirstFragmentIfKnown(parameters, proteinConformers, molecule, averageEnergy, "conformers", errorFF);
            }
        }

    // If we want to grow molecules, here is the core of the program.
    else {
       for (int x=0 ; x < parameters.numberOutputs ; x++) {
        // Define some parameters.
        int counterIterations=0;
        double averageEnergy=0;
        int successfulGrowth=1;
        int errorFF=0;
        Molecule moleculeRot[parameters.rotamersNumber];
        Molecule moleculeConf[parameters.conformersNumber];

        // Look for a value of "x" which has not been used yet.
        int xCanBeUsed=0;
        while (xCanBeUsed==0) {
            // Counters for the output file.
            ostringstream counterOutputFlow;
            counterOutputFlow << x;
            string outputNumber = counterOutputFlow.str();
            // Full output name.
            string outputFileName = parameters.outputName + "_" + outputNumber + "_0.xyz";

            // Get the size of the output file.
            int sizeOutputFile=SizeFile(outputFileName);
            // If the file is not empty, increment x; if it is empty write in it to block the "x" value.
            if   (sizeOutputFile!=0) { x++; }
            else { 
                   // Open the output file.
                   ofstream outputFile(outputFileName.c_str(), ios::out);
                   if(!outputFile) {
                       cerr << "ERROR: I cannot write in the output file " << outputFileName << "." << endl;
                       exit (-1);
                       }
                   outputFile << "  " << endl;
                   outputFile.close();
                   xCanBeUsed=1;
                 }
            }
        cout << endl << "==== Ligand " << x << " ====" << endl;

        // We will need two randoms points for the position and orientation of the first ligands. We need to define them here because the same points will be used for both the rotamer and the conformer search.
        double randomAxe[3]   = {0, 0, 0};
        double randomPoint[3] = {0, 0, 0};

        // Perform a first growth if there are some rotamers.
        if (parameters.rotamersNumber != 0) {
            if (parameters.verbose>=3) { cout << "***** Start rotamers *****" << endl; }
            // Either add the first fragment randomly in the active site or put a given one.
            if (parameters.verbose>=3) { cout << "\t\t--- Try to add fragment 1." << endl; }
            if      (parameters.mode=="SEED")   { FirstFragmentIfKnown(parameters, proteinRotamers, moleculeRot, averageEnergy, "rotamers", errorFF);                                                               }
            else if (parameters.mode=="DENOVO") { FirstFragment(parameters, proteinRotamers, moleculeRot, averageEnergy, "rotamers", fragment, probaFirstFrag, randomAxe, randomPoint, ligandDescription, errorFF); }
            if (parameters.verbose>=5) { CheckGrowth(moleculeRot[0]); }

            // Then add new fragments, as long as no threshold is reached for the growth, there are available hydrogens, and we have not tried too much iterations.
            counterIterations=0;
            while ( (moleculeRot[0].numberResidues<parameters.maxFragment) && (moleculeRot[0].obMol.NumHvyAtoms()<parameters.maxAtoms) && (moleculeRot[0].obMol.GetMolWt()<parameters.maxMW) && 
                (moleculeRot[0].growingSites.size()!=0) && (counterIterations<parameters.maxIterations) && (errorFF==0)) {
                    if (parameters.verbose>=3) { cout << "\t\t--- Try to add fragment " << moleculeRot[0].numberResidues+1 << "." << endl; }
                    int currentResiduesNumber=moleculeRot[0].numberResidues;
                    AddFragment(parameters, proteinRotamers, moleculeRot, averageEnergy, "rotamers", fragment, probaFirstFrag, probaTransition, ligandDescription, errorFF);
                    counterIterations++;
                    if (parameters.verbose>=3) { cout << "\t\t" << moleculeRot[0].obMol.NumAtoms() << " atoms, " << moleculeRot[0].growingSites.size() << " growing sites, " << moleculeRot[0].numberResidues << " residues." << endl; }
                    if (moleculeRot[0].numberResidues>currentResiduesNumber) {
                    counterIterations=0;
                    if (parameters.verbose>=3) { cout << "\t\t+++ Fragment " << moleculeRot[0].numberResidues << " successfully added." << endl; }
                    if (parameters.verbose>=5) { CheckGrowth(moleculeRot[0]); }
                    }
                }

            // Define the growth as successful if the numbers of heavy atoms and fragments are above given limits and the average energy below another.
            if ( (moleculeRot[0].obMol.NumHvyAtoms()>=parameters.minAtoms) && (moleculeRot[0].numberResidues>=parameters.minFragments) && (averageEnergy<=parameters.minEnergy) && errorFF==0 ) {
                successfulGrowth=1;
                // If there are no conformers, write the output and save the molecules because the loop will stop here.
                if ( parameters.conformersNumber==0 && errorFF==0 ) {
                    SaveOutput(parameters, proteinRotamers, moleculeRot, ligandDescription, averageEnergy, x, threeMerList, "rotamers", errorFF);
                    cout << "\tLigand " << x << " saved with a mass of " << moleculeRot[0].obMol.GetMolWt() << " g/mol" << endl;
                    }
                }
            else    { successfulGrowth = 0; }
            }

        // Then perform a growth with conformers.
        if ( (parameters.conformersNumber != 0) && (successfulGrowth == 1) ) {
            if (parameters.verbose>=3) { cout << "***** Start conformers *****" << endl; }
            // Either add the first fragment randomly in the active site or put a given one.
            if (parameters.verbose>=3) { cout << "\t\t--- Try to add fragment 1." << endl; }
            if      (parameters.mode=="SEED")   { FirstFragmentIfKnown(parameters, proteinConformers, moleculeConf, averageEnergy, "conformers", errorFF);                                                               }
            else if (parameters.mode=="DENOVO") { FirstFragment(parameters, proteinConformers, moleculeConf, averageEnergy, "conformers", fragment, probaFirstFrag, randomAxe, randomPoint, ligandDescription, errorFF); }
            if (parameters.verbose>=5) { CheckGrowth(moleculeConf[0]); }

            // Then add new fragments, as long as no threshold is reached for the growth, there are available hydrogens, and we have not tried too much iterations.
            counterIterations=0;
            while ( (moleculeConf[0].numberResidues<parameters.maxFragment) && (moleculeConf[0].obMol.NumHvyAtoms()<parameters.maxAtoms) && (moleculeConf[0].obMol.GetMolWt()<parameters.maxMW) &&
                (moleculeConf[0].growingSites.size()!=0) && (counterIterations<parameters.maxIterations) && (errorFF==0)) {
                    if (parameters.verbose>=3) { cout << "\t\t--- Try to add fragment " << moleculeConf[0].numberResidues+1 << "." << endl; }
                    int currentResiduesNumber=moleculeConf[0].numberResidues;
                    AddFragment(parameters, proteinConformers, moleculeConf, averageEnergy, "conformers", fragment, probaFirstFrag, probaTransition, ligandDescription, errorFF);
                    counterIterations++;
                    if (parameters.verbose>=3) { cout << "\t\t" << moleculeConf[0].obMol.NumAtoms() << " atoms, " << moleculeConf[0].growingSites.size() << " growing sites, " << moleculeConf[0].numberResidues << " residues." << endl; }
                    if (moleculeConf[0].numberResidues>currentResiduesNumber) {
                    counterIterations=0;
                    if (parameters.verbose>=3) { cout << "\t\t+++ Fragment " << moleculeConf[0].numberResidues << " successfully added." << endl; }
                    if (parameters.verbose>=5) { CheckGrowth(moleculeConf[0]); }
                    }
                }

            // Save the molecule if the numbers of heavy atoms and fragments are above given limits and the average energy below another.
            if ( (moleculeConf[0].obMol.NumHvyAtoms()>=parameters.minAtoms) && (moleculeConf[0].numberResidues>=parameters.minFragments) && (averageEnergy<=parameters.minEnergy) && errorFF==0 ) {
                SaveOutput(parameters, proteinConformers, moleculeConf, ligandDescription, averageEnergy, x, threeMerList, "conformers", errorFF);
                cout << "\tLigand " << x << " saved with a mass of " << moleculeConf[0].obMol.GetMolWt() << " g/mol" << endl;
                }
            else    { successfulGrowth = 0; }
            }

        // If the growth was not successfull with this "x" value, delete the file.
        if (successfulGrowth==0) {
            // Counters for the output file.
            ostringstream counterOutputFlow;
            counterOutputFlow << x;
            string outputNumber = counterOutputFlow.str();
            // Full output name.
            string outputFileName = parameters.outputName + "_" + outputNumber + "_0.xyz";
            // Remove file.
            remove(outputFileName.c_str());
            }
        } // End "for (int x=0 ; x < parameters.numberOutputs ; x++)"
    } // End "else"

    cout << endl << "######################################################################" << endl;
    cout << endl;
    cout << "**********************************************************************" << endl;
    cout << "**************              Thank you for               **************" << endl;
    cout << "**************            using OpenGrowth.             **************" << endl;
    cout << "**********************************************************************" << endl;
    cout << endl;
    return 0;
}

