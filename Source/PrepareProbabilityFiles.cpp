#include "OpenGrowth.h"

int SizeFile(string const fileName);

// This function prepares the probability files by reading files and storing the numbers.
void PrepareProbabilityFiles(Parameters & parameters, double probaFirstFrag[], double probaTransition[][MAX_FRAGMENTS])
{
    // Check that the sizes are consistent.
    if ( SizeFile(parameters.probaFirstFrag)!=parameters.fragmentListSize ) {
        cerr << "ERROR: The fragment list file and the probability file for the first fragment don't have the same size." << endl;
        exit (-1);
        }
    if ( SizeFile(parameters.probaTransition)!=parameters.fragmentListSize ) {
        cerr << "ERROR: The fragment list file and the probability file for the transition to rings don't have the same size." << endl;
        exit (-1);
        }

    // Initialize the arrays with 0.0.
    for (int j=0; j<parameters.fragmentListSize ; j++ ) {
        probaFirstFrag[j]=0.0;
        for (int i=0; i<MAX_FRAGMENTS ; i++ ) {
            probaTransition[j][i]=0.0;
            }
        }

    // Open the first fragment probability file, check if the file has been opened.
    ifstream probaFirstFragFile(parameters.probaFirstFrag.c_str(), ios::in);
    if(!probaFirstFragFile) {
        cerr << "ERROR: The file " << parameters.probaFirstFrag << " is missing, or I can't read it." << endl;
        exit (-1);   
        }
    string line;
    int index=0;
    while(getline(probaFirstFragFile, line)) {
        istringstream input(line);
        double probaValue=0.0;
        input >> probaValue;
        probaFirstFrag[index] = probaValue;
        index++;
        }
    probaFirstFragFile.close();

    // Open the transition probability file, check if the file has been opened. Do it only if GROWTH_MODE=FOG.
    if (parameters.growthMode=="FOG") {
        ifstream probaTransitionFile(parameters.probaTransition.c_str(), ios::in);
        if(!probaTransitionFile) {
            cerr << "ERROR: The file " << parameters.probaTransition << " is missing, or I can't read it." << endl;
            exit (-1);   
            }
        int j=0;
        while(getline(probaTransitionFile, line)) {
            istringstream input(line);
            for (int i=0 ; i < parameters.fragmentListSize ; i++) {
                double probaValue=0.0;
                input >> probaValue;
                probaTransition[j][i] = probaValue;
                }
            j++;
            }
        probaTransitionFile.close();
        }

    // For each fragment, check that some dimers exist: the sum of all the probabilities must be 1 and not 0.0. If some dimers doesn't exist, avoid using it as first fragments.
    if (parameters.growthMode!="RANDOM") {
        for (int j=0 ; j < parameters.fragmentListSize ; j++) {
            double Sum=0.0;
            for (int i=0 ; i < parameters.fragmentListSize ; i++) {
                Sum += probaTransition[j][i];
                }
            if (Sum==0.0) {
                probaFirstFrag[j]=0.0;
                parameters.forbiddenFragments.push_back(j+1);
                }
            }
        }

    cout << endl << "*************** Reading of probability files successful **************" << endl;
}

