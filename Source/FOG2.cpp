#include "OpenGrowth.h"

void  Parse(int nbarg, char * argv[], Parameters & parameters);
void  PrepareFragments(Parameters & parameters, Molecule fragment[]);
void  Prepare3Mer(Parameters const & parameters, string threeMerList[]);
void  PrepareProbabilityFiles(Parameters & parameters, double probaFirstFrag[], double probaTransition[][MAX_FRAGMENTS]);
int   Random(Parameters & parameters, double const probaFile[]);
float Optimization(Parameters const & parameters, Molecule & molecule, int & errorFF, int const & optParameter=1);
void  AddFragment(Parameters & parameters, Molecule const fragment[], double const probaTransition[][MAX_FRAGMENTS], double const probaFirstFrag[], Molecule & molecule, int & errorFF);
void  SaveOutput(Parameters const & parameters, Molecule & molecule, int const & counterOutput, string const threeMerList[]);
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
    ifstream listFile(fileName.c_str(), ios::in);
    if(!listFile) {
        cerr << "ERROR: The file " << fileName << " is missing, or I can't read it." << endl;
        exit (-1);
        }
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
    cout << "**************************************************************" << endl;
    cout << "**********                   FOG                    **********" << endl;
    cout << "**********                   2.0                    **********" << endl;
    cout << "**************************************************************" << endl;
    cout << "**********          Author: Nicolas Chéron          **********" << endl;
    cout << "**********         (FOG1: Peter Kutchukian)         **********" << endl;
    cout << "**********         Eugene Shakhnovich Group         **********" << endl;
    cout << "**********          Harvard CCB, 2015-2017          **********" << endl;
    cout << "**********           License: GNU GPL v3            **********" << endl;
    cout << "**************************************************************" << endl;
    cout << endl;

    // Initialize random numbers.
    srand(time(NULL));

    // Parse the input file.
    Parameters parameters;
    Parse(nbarg, argv, parameters);

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

    // Prepare fragments.
    Molecule fragment[MAX_FRAGMENTS];            // Array of fragments
    PrepareFragments(parameters, fragment);      // Prepare the fragments (store the adress, name, ring information, growth mode of atoms)

    // Prepare the forbidden 3-mer list.
    if (parameters.threeMerFile!="") { parameters.threeMerSize=SizeFile(parameters.threeMerFile); }
    else                            { parameters.threeMerSize=0; }
    string threeMerList[parameters.threeMerSize];
    Prepare3Mer(parameters, threeMerList);

    // Store information from the probability libraries.
    double probaFirstFrag[parameters.fragmentListSize];
    double probaTransition[parameters.fragmentListSize][MAX_FRAGMENTS];
    if (parameters.growthMode=="FOG" || parameters.growthMode=="BIASED") { PrepareProbabilityFiles(parameters, probaFirstFrag, probaTransition); }

    // Prepare a summary file.
    string outputSummary = parameters.outputName + "_Summary.txt";
    ofstream outputSum(outputSummary.c_str(), ios::out);
    if(!outputSum) {
        cerr << "ERROR: I cannot write in the output file " << outputSummary << "." << endl;
        exit (-1);
        }
    outputSum << "#LigandName" ;
    for (int unsigned i=0 ; i < parameters.outputName.size() ; i++) { outputSum << " "; }
    outputSum << "Fragments     Heavy atoms     Molecular weight     " ;
    if (parameters.threeMerFile!="") { outputSum << "3Mer     " ; }
    outputSum << "\tSMILES" << endl;
    outputSum.close();

    cout << endl << "######################################################################" << endl;
    if (parameters.forbiddenFragments.size()!=0) {
        cout << "The following fragments will not be used: " ;
        for (unsigned int j=0; j<parameters.forbiddenFragments.size() ; j++ ) {
            cout << parameters.forbiddenFragments[j] << " " ;
            }
        cout << endl << "######################################################################" << endl;
        }

    // Here is the core of the program.
    for (int x=0 ; x < parameters.numberOutputs ; x++) {
        // Define some parameters.
        int errorFF=0;
        Molecule molecule;
        cout << endl << "==== Molecule " << x << " ====" << endl;

        // Prepare first fragment. Choose the first fragment, either randomly or according to the FOG probabilities.
        int indexFragment=0;
        if      (parameters.growthMode=="RANDOM") { indexFragment = (parameters.random() % parameters.fragmentListSize)+1;  }
        else if (parameters.growthMode=="BIASED") { indexFragment = Random(parameters, probaFirstFrag)+1;                   }
        else if (parameters.growthMode=="FOG")    { indexFragment = Random(parameters, probaFirstFrag)+1;                   }
        molecule=fragment[indexFragment-1];

        if (parameters.verbose>=1) { cout << "Fragment " << molecule.numberResidues << ": #" << indexFragment << ", " << molecule.obMol.GetTitle() << endl; }
        if (parameters.verbose>=3) { cout << "\t\t" << molecule.obMol.NumAtoms() << " atoms, " << molecule.growingSites.size() << " growing sites, " << molecule.numberResidues << " residues." << endl; }
        Optimization(parameters, molecule, errorFF, 1);
        // Check growth.
        if (parameters.verbose>=5) { CheckGrowth(molecule); }

        // Then add new fragments as long as no treshlod is reached for the growth, there are available hydrogens, no error was found and we have not tried too many times.
        int counterIterations=0;
        while ( (molecule.numberResidues<parameters.maxFragment) && (molecule.obMol.NumHvyAtoms()<parameters.maxAtoms) && (molecule.obMol.GetMolWt()<parameters.maxMW) &&
            (molecule.growingSites.size() != 0) && (errorFF==0) && counterIterations<parameters.maxIterations ) {
                if (parameters.verbose>=3) { cout << "\t\t--- Try to add fragment " << molecule.numberResidues+1 << "." << endl; }
                int currentResiduesNumber=molecule.numberResidues;
                AddFragment(parameters, fragment, probaTransition, probaFirstFrag, molecule, errorFF);
                counterIterations++;
                if (parameters.verbose>=3) { cout << "\t\t" << molecule.obMol.NumAtoms() << " atoms, " << molecule.growingSites.size() << " growing sites, " << molecule.numberResidues << " residues." << endl; }
                if (parameters.verbose>=3 && currentResiduesNumber<molecule.numberResidues) {
                    cout << "\t\t+++ Fragment " << molecule.numberResidues << " successfully added." << endl;
                    counterIterations=0;
                    // Check growth.
                    if (parameters.verbose>=5) { CheckGrowth(molecule); }
                    }
                }

        // Save the molecule.
        if (errorFF==0) {
            Optimization(parameters, molecule, errorFF, 2);
            if ( (molecule.obMol.NumHvyAtoms()>=parameters.minAtoms) && (molecule.numberResidues>=parameters.minFragments) ) {
                SaveOutput(parameters, molecule, x, threeMerList);
                cout << "\t== Molecule " << x << " saved with a mass of " << molecule.obMol.GetMolWt() << " g/mol" << endl;
                }
            }

        molecule.obMol.Clear();
        molecule.growth.clear();
        molecule.fragmentIndex.clear();
        molecule.fragmentNeighbours.clear();
        molecule.growingSites.clear();
        molecule.numberResidues=0;
        molecule.stericClashes=0;
        }

    cout << endl ;
    cout << "##############################################################" << endl;
    cout << endl;
    cout << "**************************************************************" << endl;
    cout << "**********                Thank you for             **********" << endl;
    cout << "**********                using FOG2.0              **********" << endl;
    cout << "**************************************************************" << endl;
    cout << endl;
    return 0;
}

