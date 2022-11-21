#include "OpenGrowth.h"

int SizeFile(string const fileName);

// This function reads the parameter file and store the information.
void Parse(int nbarg, char * argv[], Parameters & parameters)
{
    // Display the help.
    if(nbarg==1) {
        cout << "The use of OpenGrowth is: ./OpenGrowth.exe InputFile.inp. The input file InputFile.inp is made of lines with \"PARAMETERNAME parametervalue\"." << endl;
        cout << "Only one parametervalue per line is read. The two entries (name and value) can be separated by spaces or tabulations." << endl;
        cout << "All parameters are mandatory, unless a default value is defined. See the documentation for more details." << endl;
        cout << endl;
        exit(0);
        }
    for(int i=1; i<nbarg; i++) {
        if (!strcmp(argv[i],"-help") || !strcmp(argv[i],"--help") || !strcmp(argv[i],"-h") || !strcmp(argv[i],"--h") || !strcmp(argv[i],"-H") || !strcmp(argv[i],"--H")) {
            cout << "The use of OpenGrowth is: ./OpenGrowth.exe InputFile.inp. The input file InputFile.inp is made of lines with \"PARAMETERNAME parametervalue\"." << endl;
            cout << "Only one parametervalue per line is read. The two entries (name and value) can be separated by spaces or tabulations." << endl;
            cout << "All parameters are mandatory, unless a default value is defined. See the documentation for more details." << endl;
            cout << endl;
            exit(0);
            }
        }

    // Open the input file.
    ifstream parameterFile(argv[1], ios::in);
    if(!parameterFile) {
        cerr << "ERROR: The file " << argv[1] << " is missing, or I can't read it." << endl;
        exit (-1);
        }
    cout << "***********   Start reading the parameter file " << argv[1] << "   ************" << endl;

    // Initialize all the parameters to 0 (or something else if 0 is one of the options) or empty strings.
    parameters.rotamersName="";
    parameters.conformersName="";
    parameters.rotamersNumber=0;
    parameters.conformersNumber=0;
    parameters.proteinRange=0;
    parameters.mode="";
    parameters.growthMode="";
    parameters.probaFirstFrag="";
    parameters.probaTransition="";
    parameters.regrowFile="";
    parameters.scoringFunction="";
    parameters.energyFile="";
    parameters.averageType="";
    parameters.ligandName="";
    parameters.bindingSite[0]=0;
    parameters.bindingSite[1]=0;
    parameters.bindingSite[2]=0;
    parameters.bindingSize=0;
    parameters.fragmentListName="";
    parameters.fragmentListSize=0;
    parameters.vdWScaleInter=0;
    parameters.vdWScaleIntra=0;
    parameters.branchingProba=10.0;
    parameters.rotationPrecision=0;
    parameters.optimizationMode=10;
    parameters.optimizationNumber=0;
    parameters.optimizationIterations=0;
    parameters.optimizationDistance=0;
    parameters.optimizationAngle=0;
    parameters.optimizationForceField="";
    parameters.optimizationSteepestDescent=99999;
    parameters.optimizationConjugateGradient=99999;
    parameters.optimizationVdwCutOff=100.0;
    parameters.optimizationElecCutOff=100.0;
    parameters.MCTemp=0;
    parameters.maxFragment=0;
    parameters.maxAtoms=0;
    parameters.maxMW=0;
    parameters.maxIterations=0;
    parameters.minFragments=0;
    parameters.minAtoms=0;
    parameters.minEnergy=0;
    parameters.outputName="";
    parameters.smilesOnly=0;
    parameters.writeDescription=0;
    parameters.threeMerFile="";
    parameters.threeMerSize=0;
    parameters.verbose=10;
    parameters.numberOutputs=0;
    parameters.forbiddenFragments.clear();

//////////////////////

    // Parse the parameterFile by looking for specific words. Look for PARAMETERNAME keywords, then store the parameter value.
    string line;
    while(getline(parameterFile, line)) {
        istringstream input(line);
        string parameterName="";
        string parameterValue="";
        input >> parameterName >> parameterValue;
        // Receptor
        if      (parameterName=="ROTAMERS") {
            parameters.rotamersName = parameterValue;
            cout << "The rotamer name is:                                          " << parameters.rotamersName << endl;
            }
        else if (parameterName=="CONFORMERS") {
            parameters.conformersName = parameterValue;
            cout << "The conformer name is:                                        " << parameters.conformersName << endl;
            }
        else if (parameterName=="ROTAMERS_NUMBER") {
            parameters.rotamersNumber=atoi(parameterValue.c_str());
            if (parameters.rotamersNumber!=0) { cout << "Number of rotamers to be used:                                " << parameters.rotamersNumber << endl; }
            }
        else if (parameterName=="CONFORMERS_NUMBER") {
            parameters.conformersNumber=atoi(parameterValue.c_str());
            cout << "Number of conformers to be used:                              " << parameters.conformersNumber << endl;
            }
        else if (parameterName=="PROTEIN_RANGE") {
            parameters.proteinRange=atof(parameterValue.c_str());
            cout << "The protein range is:                                         " << parameters.proteinRange << endl;
            }
        // Growing
        else if (parameterName=="MODE") {
            parameters.mode = parameterValue;
            cout << "The program mode is:                                          " << parameters.mode << endl;
            }
        else if (parameterName=="GROWTH_MODE") {
            parameters.growthMode = parameterValue;
            cout << "The growth mode is:                                           " << parameters.growthMode << endl;
            }
        else if (parameterName=="PROBA_FIRSTFRAG") {
            parameters.probaFirstFrag = parameterValue;
            cout << "The first fragment probability file is:                       " << parameters.probaFirstFrag << endl;
            }
        else if (parameterName=="PROBA_TRANSITION") {
            parameters.probaTransition = parameterValue;
            cout << "The transition probability file is:                           " << parameters.probaTransition << endl;
            }
        else if (parameterName=="REGROW_FILE") {
            parameters.regrowFile = parameterValue;
            cout << "The regrow file is:                                           " << parameters.regrowFile << endl;
            }
        // Scoring
        else if (parameterName=="SCORING_FUNCTION") {
            parameters.scoringFunction = parameterValue;
            cout << "The scoring function is:                                      " << parameters.scoringFunction << endl;
            }
        else if (parameterName=="ENERGY_FILE") {
            parameters.energyFile = parameterValue;
            cout << "The energetic parameters will be read from the file:          " << parameters.energyFile << endl;
            }
        else if (parameterName=="AVERAGE_TYPE") {
            parameters.averageType = parameterValue;
            cout << "The type of average is:                                       " << parameters.averageType << endl;
            }
        else if (parameterName=="LIGAND") {
            parameters.ligandName = parameterValue;
            cout << "The ligand name is:                                           " << parameters.ligandName << endl;
            }
        // Active site
        else if (parameterName=="BINDING_SITE_X") {
            parameters.bindingSite[0]=atof(parameterValue.c_str());
            cout << "The X binding site is:                                        " << parameters.bindingSite[0] << endl;
            }
        else if (parameterName=="BINDING_SITE_Y") {
            parameters.bindingSite[1]=atof(parameterValue.c_str());
            cout << "The Y binding site is:                                        " << parameters.bindingSite[1] << endl;
            }
        else if (parameterName=="BINDING_SITE_Z") {
            parameters.bindingSite[2]=atof(parameterValue.c_str());
            cout << "The Z binding site is:                                        " << parameters.bindingSite[2] << endl;
            }
        else if (parameterName=="BINDINGBOX_SIZE") {
            parameters.bindingSize=atof(parameterValue.c_str());
            cout << "The binding box size is:                                      " << parameters.bindingSize << endl;
            }
        else if (parameterName=="FRAGMENT_LIST") {
            parameters.fragmentListName = parameterValue;
            cout << "The fragment list file is:                                    " << parameters.fragmentListName << endl;
            }
        // Optimization
        else if (parameterName=="VDW_SCALE_INTER") {
            parameters.vdWScaleInter=atof(parameterValue.c_str());
            cout << "The vdW scale parameter for inter steric clashes is:          " << parameters.vdWScaleInter << endl;
            }
        else if (parameterName=="VDW_SCALE_INTRA") {
            parameters.vdWScaleIntra=atof(parameterValue.c_str());
            cout << "The vdW scale parameter for intra steric clashes is:          " << parameters.vdWScaleIntra << endl;
            }
        else if (parameterName=="BRANCHING_PROBA") {
            parameters.branchingProba=atof(parameterValue.c_str());
            cout << "The branching probability is:                                 " << parameters.branchingProba << endl;
            }
        else if (parameterName=="ROTATION_PRECISION") {
            parameters.rotationPrecision=atoi(parameterValue.c_str());
            cout << "The rotation precision is:                                    " << parameters.rotationPrecision << endl;
            }
        else if (parameterName=="OPTIMIZATION_MODE") {
            parameters.optimizationMode=atoi(parameterValue.c_str());
            cout << "The mode for optimization is:                                 " << parameters.optimizationMode << endl;
            }
        else if (parameterName=="OPTIMIZATION_NUMBER") {
            parameters.optimizationNumber=atoi(parameterValue.c_str());
            cout << "The number of optimization iterations is:                     " << parameters.optimizationNumber << endl;
            }
        else if (parameterName=="OPTIMIZATION_ITERATIONS") {
            parameters.optimizationIterations=atoi(parameterValue.c_str());
            cout << "The number of optimization per iteration is:                  " << parameters.optimizationIterations << endl;
            }
        else if (parameterName=="OPTIMIZATION_DISTANCE") {
            parameters.optimizationDistance=atof(parameterValue.c_str());
            cout << "The step for translation is:                                  " << parameters.optimizationDistance << endl;
            }
        else if (parameterName=="OPTIMIZATION_ANGLE") {
            parameters.optimizationAngle=atof(parameterValue.c_str());
            cout << "The angle for rotation is:                                    " << parameters.optimizationAngle << endl;
            }
        // Force field for optimization
        else if (parameterName=="OPTIMIZATION_FORCEFIELD") {
            parameters.optimizationForceField = parameterValue;
            cout << "The force field used for optimization is:                     " << parameters.optimizationForceField << endl;
            }
        else if (parameterName=="OPTIMIZATION_STEEPDESC") {
            parameters.optimizationSteepestDescent=atoi(parameterValue.c_str());
            cout << "The number of steepest descent steps performed is:            " << parameters.optimizationSteepestDescent << endl;
            }
        else if (parameterName=="OPTIMIZATION_CONJGRAD") {
            parameters.optimizationConjugateGradient=atoi(parameterValue.c_str());
            cout << "The number of conjugate gradient steps performed is:          " << parameters.optimizationConjugateGradient << endl;
            }
        else if (parameterName=="OPTIMIZATION_VDWCUTOFF") {
            parameters.optimizationVdwCutOff=atof(parameterValue.c_str());
            cout << "The cut-off for vdW interactions is:                          " << parameters.optimizationVdwCutOff << endl;
            }
        else if (parameterName=="OPTIMIZATION_ELECCUTOFF") {
            parameters.optimizationElecCutOff=atof(parameterValue.c_str());
            cout << "The cut-off for electrostatic interactions is:                " << parameters.optimizationElecCutOff << endl;
            }
        else if (parameterName=="MC_TEMP") {
            parameters.MCTemp=atof(parameterValue.c_str());
            cout << "The Monte Carlo temperature is:                               " << parameters.MCTemp << endl;
            }
        // Threshold
        else if (parameterName=="MAX_FRAGMENTS") {
            parameters.maxFragment=atoi(parameterValue.c_str());
            cout << "The maximum number of fragments is:                           " << parameters.maxFragment << endl;
            }
        else if (parameterName=="MAX_ATOMS") {
            parameters.maxAtoms=atoi(parameterValue.c_str());
            cout << "The maximum number of heavy atoms is:                         " << parameters.maxAtoms << endl;
            }
        else if (parameterName=="MAX_MW") {
            parameters.maxMW=atoi(parameterValue.c_str());
            cout << "The growth will stop at a molecular weight of:                " << parameters.maxMW << " g/mol" << endl;
            }
        else if (parameterName=="MAX_ITERATIONS") {
            parameters.maxIterations=atoi(parameterValue.c_str());
            cout << "The maximum number of iterations is:                          " << parameters.maxIterations << endl;
            }
        // Miscellaneous
        else if (parameterName=="MIN_FRAGMENTS") {
            parameters.minFragments=atoi(parameterValue.c_str());
            cout << "The minimum number of fragments is:                           " << parameters.minFragments << endl;
            }
        else if (parameterName=="MIN_ATOMS") {
            parameters.minAtoms=atoi(parameterValue.c_str());
            cout << "The minimum number of heavy atoms is:                         " << parameters.minAtoms << endl;
            }
        else if (parameterName=="MIN_ENERGY") {
            parameters.minEnergy=atof(parameterValue.c_str());
            cout << "The minimum energy is (we keep struct. with a lower score):   " << parameters.minEnergy << endl;
            }
        else if (parameterName=="OUTPUT") {
            parameters.outputName = parameterValue;
            cout << "The output files will be called:                              \"" << parameters.outputName << "\"" << endl;
            }
        else if (parameterName=="SMILESONLY") {
            parameters.smilesOnly = atoi(parameterValue.c_str());
            cout << "The parameter for writing only SMILES string is:              " << parameters.smilesOnly << endl;
            }
        else if (parameterName=="WRITEDESCRIPTION") {
            parameters.writeDescription = atoi(parameterValue.c_str());
            cout << "The parameter for writing the ligand description is:          " << parameters.writeDescription << endl;
            }
        else if (parameterName=="3MERSCREEN") {
            parameters.threeMerFile = parameterValue;
            cout << "The 3Mer-Screen file is:                                      " << parameters.threeMerFile << endl;
            }
        else if (parameterName=="VERBOSE") {
            parameters.verbose=atoi(parameterValue.c_str());
            cout << "The verbose level is:                                         " << parameters.verbose << endl;
            }
        else if (parameterName=="NUMBER_OUTPUT") {
            parameters.numberOutputs=atoi(parameterValue.c_str());
            cout << "We will grow:                                                 " << parameters.numberOutputs << " ligand(s)" << endl;
            }
        }

//////////////////////

    // Check that the parameters exist. If a parameter is missing, assign the default value if there is one or stop the program otherwise.
    // Receptor
    if (parameters.rotamersName=="" && parameters.conformersName=="") {
        cerr << "ERROR: You must specify either a rotamer file name (ROTAMERS) or a conformer file name (CONFORMERS)." << endl;
        exit(-1);
        }
    if (!parameters.rotamersNumber && !parameters.conformersNumber) {
        cerr << "ERROR: You must specify the number of rotamers (ROTAMERS_NUMBER) and/or conformers (CONFORMERS_NUMBER) that will be used." << endl;
        exit(-1);
        }
    if (!parameters.rotamersNumber) {
        cout << "No parameter found for the number of rotamers (ROTAMERS_NUMBER). We will not use them." << endl;
        }
    if (!parameters.conformersNumber) {
        cout << "No parameter found for the number of conformers (CONFORMERS_NUMBER). We will not use them." << endl;
        }
    if (parameters.rotamersName=="") {
        cout << "No rotamers file name was found." << endl;
        parameters.rotamersNumber=0;
        }
    if (parameters.conformersName=="") {
        cout << "No conformers file name was found." << endl;
        parameters.conformersNumber=0;
        }
    // Scoring
    if (!parameters.proteinRange) {
        cout << "No parameter found for the protein range (PROTEIN_RANGE). Default value of 40 Angstrom will be used." << endl;
        parameters.proteinRange=40;
        }
    if (parameters.scoringFunction=="") {
        cout << "No parameter found for the scoring function (SCORING_FUNCTION). Default value SMOG2001 will be used." << endl;
        parameters.scoringFunction="SMOG2001";
        }
    if (parameters.scoringFunction!="SMOG2001" && parameters.scoringFunction!="SMOG2016") {                      // And later the other possible scoring functions.
        cerr << "ERROR: The only possible options for SCORING_FUNCTION are SMOG2001 or SMOG2016." << endl;
        exit(-1);
        }
    if (parameters.scoringFunction=="SMOG2016") {
        if (SizeFile("VDWParameters.txt")==0) {
            cout << "WARNING: The file VDWParameters.txt is missing whereas it is needed by SMOG2016. The SMOG2001 scoring function will be used." << endl;
            parameters.scoringFunction="SMOG2001";
            }
        if (SizeFile("AmberAtomTypes.txt")==0) {
            cout << "WARNING: The file AmberAtomTypes.txt is missing whereas it is needed by SMOG2016. The SMOG2001 scoring function will be used." << endl;
            parameters.scoringFunction="SMOG2001";
            }
        if (parameters.energyFile=="") {
            cout << "WARNING: No file was defined as \"ENERGY_FILE\" whereas one is needed by SMOG2016. The SMOG2001 scoring function will be used." << endl;
            parameters.scoringFunction="SMOG2001";
            }
        }

    // Active site
    if (!parameters.bindingSite[0] || !parameters.bindingSite[1] || !parameters.bindingSite[2]) {
        cerr << "ERROR: You must specify a correct binding site (BINDING_SITE_X / BINDING_SITE_Y / BINDING_SITE_Z). If one of the value is exactly 0, change it to 0.01: it will not change the science but the program will work." << endl;
        exit(-1);
        }
    if (!parameters.bindingSize) {
        cout << "No parameter found for the binding box size (BINDINGBOX_SIZE). Default value of 1.0 Angstrom will be used." << endl;
        parameters.bindingSize=1.0;
        }
    // Growing
    if (parameters.mode=="") {
        cout << "No parameter found for the program mode (MODE). Default value \"DENOVO\" will be used." << endl;
        parameters.mode="DENOVO";
        }
    if (parameters.mode!="DENOVO" && parameters.mode!="SEED" && parameters.mode!="ENERGY") {
        cerr << "ERROR: The only possible options for MODE are DENOVO, SEED, ENERGY. If nothing is selected, the default value \"DENOVO\" will be used." << endl;
        exit(-1);
        }
    if (parameters.growthMode=="") {
        cout << "No parameter found for the growth mode (GROWTH_MODE). Default value \"RANDOM\" will be used." << endl;
        parameters.growthMode="RANDOM";
        }
    if (parameters.growthMode!="RANDOM" && parameters.growthMode!="BIASED" && parameters.growthMode!="FOG" && parameters.growthMode!="REGROW") {
        cerr << "ERROR: The only possible options for GROWTH_MODE are RANDOM, BIASED, FOG, REGROW. If nothing is selected, the default value \"RANDOM\" will be used." << endl;
        exit(-1);
        }
    if (parameters.growthMode=="BIASED" && parameters.probaFirstFrag=="") {
        cout << "WARNING: No parameter found for the first fragments probabilities (PROBA_FIRSTFRAG). We will switch to a random growth." << endl;
        parameters.growthMode="RANDOM";
        }
    if (parameters.growthMode=="BIASED" && parameters.probaTransition=="") {
        parameters.probaTransition=parameters.probaFirstFrag;
        }
    if (parameters.growthMode=="FOG" && (parameters.probaFirstFrag=="" || parameters.probaTransition=="")) {
        cout << "WARNING: No parameter found for either the first fragments probabilities (PROBA_FIRSTFRAG) or the probabilities of transition (PROBA_TRANSITION). We will switch to a RANDOM growth." << endl;
        parameters.growthMode="RANDOM";
        }
    if (parameters.ligandName=="") {
        if (parameters.mode=="SEED" || parameters.mode=="ENERGY") {
            cerr << "ERROR: When MODE is SEED or ENERGY, you must specify a ligand file name (LIGAND)." << endl;
            exit(-1);
            }
        }
    if (parameters.growthMode=="REGROW" && parameters.regrowFile=="") {
        cerr << "ERROR: No parameter found for the Regrow file." << endl;
        exit(-1);
        }
    if (parameters.branchingProba==10.0) {
        cout << "No parameter found for branching probability (BRANCHING_PROBA). Default value of 0.5 will be used." << endl;
        parameters.branchingProba=0.5;
        }
    if (!parameters.MCTemp) {
        cout << "No parameter found for the Monte Carlo temperature (MC_TEMP). Default value of 1.3 will be used." << endl;
        parameters.MCTemp=1.3;
        }
    if (parameters.fragmentListName=="") {
        cerr << "ERROR: You must specify a fragment list name (FRAGMENT_LIST)." << endl;
        exit(-1);
        }
    // Optimization
    if (!parameters.vdWScaleInter) {
         cout << "No parameter found for the intermolecular vdW steric clashes parameter (VDW_SCALE_INTER). Default value of 0.70 will be used." << endl;
         parameters.vdWScaleInter=0.70;
        }
    if (!parameters.vdWScaleIntra) {
         cout << "No parameter found for the intramolecular vdW steric clashes parameter (VDW_SCALE_INTRA). Default value of 0.70 will be used." << endl;
         parameters.vdWScaleInter=0.70;
        }
	if (!parameters.rotationPrecision) {
        cout << "No parameter found for the rotation precision (ROTATION_PRECISION). Default value of 24 will be used (i.e. every 15°)." << endl;
        parameters.rotationPrecision=24;
        }
    if (parameters.optimizationMode>=100) {
        cout << "No correct parameter found for the optimization mode (OPTIMIZATION_MODE). Default value of 12 will be used." << endl;
        parameters.optimizationMode=12;
        }
    if (!parameters.optimizationNumber) {
        cout << "No parameter found for the optimization number (OPTIMIZATION_NUMBER). Default value of 30 will be used." << endl;
        parameters.optimizationNumber=30;
        }
    if (!parameters.optimizationIterations) {
        cout << "No parameter found for the iteration number (OPTIMIZATION_ITERATIONS). Default value of 10 will be used." << endl;
        parameters.optimizationIterations=10;
        }
    if (!parameters.optimizationDistance) {
        cout << "No parameter found for the translation step (OPTIMIZATION_DISTANCE). Default value of 0.05 Angstrom will be used." << endl;
        parameters.optimizationDistance=0.05;
        }
    if (!parameters.optimizationAngle) {
        cout << "No parameter found for the rotation step (OPTIMIZATION_ANGLE). Default value of 1 degree will be used." << endl;
        parameters.optimizationAngle=1;
        }
    // Force field for optimization
    if (parameters.optimizationForceField=="") {
        cout << "No parameter found for the force field for the optimization (OPTIMIZATION_FORCEFIELD). MMFF94 will be used." << endl;
        parameters.optimizationForceField="MMFF94";
        }
    if (parameters.optimizationForceField!="MMFF94" && parameters.optimizationForceField!="UFF" && parameters.optimizationForceField!="Ghemical" && parameters.optimizationForceField!="GAFF") {
        cerr << "ERROR: The only possible options for OPTIMIZATION_FORCEFIELD are MMFF94, UFF, Ghemical, GAFF. If nothing is selected, the default value \"MMFF94\" will be used." << endl;
        exit(-1);
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
    // Threshold
    if (!parameters.maxFragment && !parameters.maxAtoms && !parameters.maxMW) {
        cerr << "ERROR: Either a maximum number of fragment (MAX_FRAGMENT), a maximum number of heavy atoms (MAX_ATOMS) or a maximum molecular weight (MAX_MW) must be provided." << endl;
        exit(-1);
        }
    if (!parameters.maxIterations) {
        cout << "No parameter found for the maximum number of iterations (MAX_ITERATIONS). Default value of 20 will be used." << endl;
        parameters.maxIterations=20;
        }
    if (!parameters.minFragments) {
        cout << "No parameter found for the minimum number of fragments (MIN_FRAGMENTS). Default value of 0 will be used." << endl;
        }
    if (!parameters.minAtoms) {
        cout << "No parameter found for the minimum number of heavy atoms (MIN_ATOMS). Default value of 5 will be used." << endl;
        parameters.minAtoms=5;
        }
    if (!parameters.minEnergy) {
        cout << "No parameter found for the minimum energy (MIN_ENERGY). Default value of 0 will be used." << endl;
        }
    // Miscellaneous
    if (parameters.outputName=="") {
        cout << "No parameter found for the output file name (OUTPUT). Default value of \"Ligand\" will be used." << endl;
        parameters.outputName = "Ligand";
        }
    if(parameters.smilesOnly!=0) {
        cout << "We will only write SMILES string." << endl;
        }
    if(parameters.writeDescription!=0) {
        cout << "We will write the ligand description." << endl;
        }
    if (parameters.threeMerFile=="") {
        cout << "No 3Mer-Screen will be performed." << endl;
        }
    if (parameters.averageType=="") {
        parameters.averageType="ARITHMETIC";
        cout << "No parameter found for the average type (AVERAGE_TYPE). We will use a arithmetic average." << endl;
        }
    if (parameters.averageType!="BOLTZMANN" && parameters.averageType!="ARITHMETIC" && parameters.averageType!="LOWESTSCORE") {
        cerr << "ERROR: The only possible options for AVERAGE_TYPE are ARITHMETIC, BOLTZMANN or LOWESTSCORE." << endl;
        exit(-1);
        }
    if (!parameters.numberOutputs) {
        cout << "No parameter found for the number of ligands you want (NUMBER_OUTPUT). Default value of 1,000,000 will be used" << endl;
        parameters.numberOutputs=1000000;
        }
    if (parameters.verbose==10) {
        cout << "No parameter found for the verbose level (VERBOSE). Default value of 2 will be used." << endl;
        parameters.verbose=2;
        }
    cout << "************** Reading of the parameter file successful **************" << endl;

//////////////////////

    // During the growth of molecules, a test on the number of fragments, the number of heavy atoms and the molecular weight is made. For this test to work, the corresponding values
    // must be different than 0. So we define some parameters in case they have not been defined before. The ideas is to put thresholds at very high values so they are never reached.
    if (!parameters.maxFragment) { parameters.maxFragment = 1000;  }
    if (!parameters.maxAtoms)    { parameters.maxAtoms    = 10000; }
    if (!parameters.maxMW)       { parameters.maxMW       = 50000; }

    // Close the file.
    parameterFile.close();
}

