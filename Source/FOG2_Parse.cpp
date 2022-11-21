#include "OpenGrowth.h"

// This function reads the parameter file and store the information.
void Parse(int nbarg, char * argv[], Parameters & parameters)
{
    // Display the help.
    if(nbarg==1) {
        cout << "The use of FOG2.0 is: ./FOG2.exe InputFile.inp. The input file InputFile.inp is made of lines with \"PARAMETERNAME parametervalue\"." << endl;
        cout << "PARAMETERNAME is case sensitive. Only one parametervalue per line is read. The two entries (name and value) can be separated by" << endl;
        cout << "spaces or tabulations. All parameters are mandatory, unless a default value is defined. See the help file for more details." << endl;
        cout << endl;
        exit(0);
        }
    for(int i=1; i<nbarg; i++) {
        if (!strcmp(argv[i],"-help") || !strcmp(argv[i],"--help") || !strcmp(argv[i],"-h") || !strcmp(argv[i],"--h") || !strcmp(argv[i],"-H") || !strcmp(argv[i],"--H")) {
            cout << "The use of FOG2.0 is: ./FOG2.exe InputFile.inp. The input file InputFile.inp is made of lines with \"PARAMETERNAME parametervalue\"." << endl;
            cout << "PARAMETERNAME is case sensitive. Only one parametervalue per line is read. The two entries (name and value) can be separated by" << endl;
            cout << "spaces or tabulations. All parameters are mandatory, unless a default value is defined. See the help file for more details." << endl;
            cout << endl;
            exit(0);
            }
        }

    // Open the file.
    ifstream parameterFile(argv[1], ios::in);
    if(!parameterFile) {
        cerr << "ERROR: The file " << argv[1] << " is missing, or I can't read it." << endl;
        exit (-1);
        }
    cout << "***********   Start reading the parameter file " << argv[1] << "   ************" << endl;

    // Initialize all the parameters to 0 (or something else if 0 is one of the options) or empty strings.
    parameters.growthMode="";
    parameters.probaFirstFrag="";
    parameters.probaTransition="";
    parameters.fragmentListName="";
    parameters.fragmentListSize=0;
    parameters.threeMerFile="";
    parameters.threeMerSize=0;
    parameters.branchingProba=10.0;
    parameters.vdWScaleIntra=0;
    parameters.rotationPrecision=0;
    parameters.optimizationForceField="";
    parameters.optimizationSteepestDescent=99999;
    parameters.optimizationConjugateGradient=99999;
    parameters.optimizationVdwCutOff=100.0;
    parameters.optimizationElecCutOff=100.0;
    parameters.maxFragment=0;
    parameters.maxAtoms=0;
    parameters.maxMW=0;
    parameters.maxIterations=0;
    parameters.outputName="";
    parameters.smilesOnly=0;
    parameters.minFragments=0;
    parameters.minAtoms=0;
    parameters.verbose=10;
    parameters.numberOutputs=0;

    // Parse the parameterFile by looking for specific words. Look for PARAMETERNAME keywords, and store them.
    string line;
    while(getline(parameterFile, line)) {
        istringstream input(line);
        string parameterName="";
        string parameterKeyword="";
        input >> parameterName >> parameterKeyword;
        if     (parameterName=="GROWTH_MODE") {
            parameters.growthMode = parameterKeyword;
            cout << "The growth mode is:                                  " << parameters.growthMode << endl;
            }
        else if (parameterName=="PROBA_FIRSTFRAG") {
            parameters.probaFirstFrag = parameterKeyword;
            cout << "The first fragment probability file is:              " << parameters.probaFirstFrag << endl;
            }
        else if (parameterName=="PROBA_TRANSITION") {
            parameters.probaTransition = parameterKeyword;
            cout << "The fragment transition probability file is:         " << parameters.probaTransition << endl;
            }
        else if (parameterName=="FRAGMENT_LIST") {
            parameters.fragmentListName = parameterKeyword;
            cout << "The fragment list file is:                           " << parameters.fragmentListName << endl;
            }
        else if (parameterName=="3MERSCREEN") {
            parameters.threeMerFile = parameterKeyword;
            cout << "The 3Mer-Screen file is:                             " << parameters.threeMerFile << endl;
            }
        else if (parameterName=="BRANCHING_PROBA") {
            parameters.branchingProba = atof(parameterKeyword.c_str());
            cout << "The branching probability is:                        " << parameters.branchingProba << endl;
            }
        else if (parameterName=="VDW_SCALE_INTRA") {
            parameters.vdWScaleIntra = atof(parameterKeyword.c_str());
            cout << "The vdW parameter for intra steric clashes is:       " << parameters.vdWScaleIntra << endl;
            }
        else if (parameterName=="ROTATION_PRECISION") {
            parameters.rotationPrecision = atoi(parameterKeyword.c_str());
            cout << "The rotation precision is:                           " << parameters.rotationPrecision << endl;
            }
        else if (parameterName=="OPTIMIZATION_FORCEFIELD") {
            parameters.optimizationForceField = parameterKeyword;
            cout << "The force field used for optimization is:            " << parameters.optimizationForceField << endl;
            }
        else if (parameterName=="OPTIMIZATION_STEEPDESC") {
            parameters.optimizationSteepestDescent = atoi(parameterKeyword.c_str());
            cout << "The number of steepest descent steps performed is:   " << parameters.optimizationSteepestDescent << endl;
            }
        else if (parameterName=="OPTIMIZATION_CONJGRAD") {
            parameters.optimizationConjugateGradient = atoi(parameterKeyword.c_str());
            cout << "The number of conjugate gradient steps performed is: " << parameters.optimizationConjugateGradient << endl;
            }
        else if (parameterName=="OPTIMIZATION_VDWCUTOFF") {
            parameters.optimizationVdwCutOff = atof(parameterKeyword.c_str());
            cout << "The cut-off for vdW interactions is:                 " << parameters.optimizationVdwCutOff << endl;
            }
        else if (parameterName=="OPTIMIZATION_ELECCUTOFF") {
            parameters.optimizationElecCutOff = atof(parameterKeyword.c_str());
            cout << "The cut-off for electrostatic interactions is:       " << parameters.optimizationElecCutOff << endl;
            }
        else if (parameterName=="MAX_FRAGMENTS") {
            parameters.maxFragment = atoi(parameterKeyword.c_str());
            cout << "The maximum number of fragments is:                  " << parameters.maxFragment << endl;
            }
        else if (parameterName=="MAX_ATOMS") {
            parameters.maxAtoms = atoi(parameterKeyword.c_str());
            cout << "The maximum number of heavy atoms is:                " << parameters.maxAtoms << endl;
            }
        else if (parameterName=="MAX_MW") {
            parameters.maxMW = atoi(parameterKeyword.c_str());
            cout << "The growth will stop at a molecular weight of:       " << parameters.maxMW << " g/mol" << endl;
            }
        else if (parameterName=="MAX_ITERATIONS") {
            parameters.maxIterations = atoi(parameterKeyword.c_str());
            cout << "The maximum number of iterations is:                 " << parameters.maxIterations << endl;
            }
        else if (parameterName=="OUTPUT") {
            parameters.outputName = parameterKeyword;
            cout << "The output files will be called:                     \"" << parameters.outputName << "\"" << endl;
            }
        else if (parameterName=="SMILESONLY") {
            parameters.smilesOnly = atoi(parameterKeyword.c_str());
            cout << "The parameter for writing only SMILES string is:     " << parameters.smilesOnly << endl;
            }
        else if (parameterName=="MIN_FRAGMENTS") {
            parameters.minFragments = atoi(parameterKeyword.c_str());
            cout << "The minimum number of fragments is:                  " << parameters.minFragments << endl;
            }
        else if (parameterName=="MIN_ATOMS") {
            parameters.minAtoms = atoi(parameterKeyword.c_str());
            cout << "The minimum number of heavy atoms is:                " << parameters.minAtoms << endl;
            }
        else if (parameterName=="VERBOSE") {
            parameters.verbose = atoi(parameterKeyword.c_str());
            cout << "The verbose level is:                                " << parameters.verbose << endl;
            }
        else if (parameterName=="NUMBER_OUTPUT") {
            parameters.numberOutputs = atoi(parameterKeyword.c_str());
            cout << "We will grow:                                        " << parameters.numberOutputs << " ligand(s)" << endl;
            }
        }

    // Check that the parameters exist. If a parameter is missing, assign the default value if there is one or stop the program otherwise.
    if (parameters.growthMode=="") {
        cout << "No parameter found for the growth mode (GROWTH_MODE). Default value \"RANDOM\" will be used." << endl;
        parameters.growthMode="RANDOM";
        }
    if (parameters.growthMode!="RANDOM" && parameters.growthMode!="BIASED" && parameters.growthMode!="FOG") {
        cerr << "ERROR: The only possible options for GROWTH_MODE are RANDOM, BIASED or FOG. If nothing is selected, the default value \"RANDOM\" will be used." << endl;
        exit(-1);
        }
    if (parameters.growthMode=="BIASED" && parameters.probaFirstFrag=="") {
        cout << "No parameter found for the first fragments probabilities (PROBA_FIRSTFRAG). We will switch to a random growth." << endl;
        parameters.growthMode="RANDOM";
        }
    if (parameters.growthMode=="BIASED" && parameters.probaTransition=="") {
        parameters.probaTransition=parameters.probaFirstFrag;
        }
    if (parameters.growthMode=="FOG" && parameters.probaFirstFrag=="" && parameters.probaTransition=="") {
        cout << "No parameter found for either the first fragments probabilities (PROBA_FIRSTFRAG) or the probabilities of transition (PROBA_TRANSITION). We will switch to a random growth." << endl;
        parameters.growthMode="RANDOM";
        }
    if (parameters.fragmentListName=="") {
        cerr << "ERROR: You must specify a fragment list name (FRAGMENT_LIST)." << endl;
        exit(-1);
        }
    if (parameters.threeMerFile=="") {
        cout << "No 3Mer-Screen will be performed." << endl;
        }
    if (parameters.branchingProba==10.0) {
        cout << "No parameter found for branching probability (BRANCHING_PROBA). Default value of 0.5 will be used." << endl;
        parameters.branchingProba=0.5;
        }
    if (!parameters.vdWScaleIntra) {
        cout << "No parameter found for the intramolecular vdW steric clashes parameter (VDW_SCALE_INTRA). Default value of 0.75 will be used." << endl;
        parameters.vdWScaleIntra=0.75;
        }
    if (!parameters.rotationPrecision) {
        cout << "No parameter found for the rotation precision (ROTATION_PRECISION). Default value of 24 will be used." << endl;
        parameters.rotationPrecision=24;
        }
    if (parameters.optimizationForceField=="") {
        cout << "No parameter found for force field name for the optimization (OPTIMIZATION_FORCEFIELD). MMFF94 will be used." << endl;
        parameters.optimizationForceField="MMFF94";
        }
    if (parameters.optimizationSteepestDescent==99999) {
        cout << "No parameter found for number of steepest descent steps (OPTIMIZATION_STEEPDESC). Default value of 100 will be used." << endl;
        parameters.optimizationSteepestDescent=100;
        }
    if (parameters.optimizationConjugateGradient==99999) {
        cout << "No parameter found for number of conjugate gradient steps (OPTIMIZATION_CONJGRAD). Default value of 0 will be used." << endl;
        parameters.optimizationConjugateGradient=0;
        }
    if (parameters.optimizationVdwCutOff==100.0) {
        cout << "No parameter found for the vdW interactions cut-off (OPTIMIZATION_VDWCUTOFF). Default value of 7.0 Angstrom will be used." << endl;
        parameters.optimizationVdwCutOff=7.0;
        }
    if (parameters.optimizationElecCutOff==100.0) {
        cout << "No parameter found for the electrostatic interactions cut-off (OPTIMIZATION_ELECCUTOFF). Default value of 15.0 Angstrom will be used." << endl;
        parameters.optimizationElecCutOff=15.0;
        }
    if (!parameters.maxFragment && !parameters.maxAtoms && !parameters.maxMW) {
        cerr << "ERROR: Either a max number of fragment (MAX_FRAGMENT), a max number of heavy atoms (MAX_ATOMS) or a max molecular weight (MAX_MW) must be provided." << endl;
        exit(-1);
        }
    if (!parameters.maxIterations) {
        cout << "No parameter found for the maximum number of iterations (MAX_ITERATIONS). Default value of 20 will be used." << endl;
        parameters.maxIterations=20;
        }
    if (parameters.outputName=="") {
        cout << "No parameter found for the output file name (OUTPUT). Default value of \"Ligand\" will be used." << endl;
        parameters.outputName="Ligand";
        }
    if(parameters.smilesOnly!=0) {
        cout << "We will only write SMILES string" << endl;
        }
    if (!parameters.minFragments) {
        cout << "No parameter found for the minimum number of fragments (MIN_FRAGMENTS). Default value of 0 will be used." << endl;
        }
    if (!parameters.minAtoms) {
        cout << "No parameter found for the minimum number of atoms (MIN_ATOMS). Default value of 5 will be used." << endl;
        parameters.minAtoms=5;
        }
    if (parameters.verbose==10) {
        cout << "No parameter found for the verbose level (VERBOSE). Default value of 2 will be used." << endl;
        parameters.verbose=2;
        }
    if (!parameters.numberOutputs) {
        cout << "No parameter found for the number of ligands you want (NUMBER_OUTPUT). Default value of 1,000,000 will be used" << endl;
        parameters.numberOutputs=1000000;
        }
    cout << "************** Reading of the parameter file successful **************" << endl;

    // During the growth of molecules, a test on the number of fragments, the number of heavy atoms and the molecular weight is made. For this test to work, the corresponding values
    // must be different than 0. So we define some parameters in case they have not been defined before. The ideas is to put tresholds at very high values so they are never reached.
    if (!parameters.maxFragment) { parameters.maxFragment = 1000;  }
    if (!parameters.maxAtoms)    { parameters.maxAtoms    = 10000; }
    if (!parameters.maxMW)       { parameters.maxMW       = 50000; }

    // Close the file.
    parameterFile.close();
}

