#include "OpenGrowth.h"

// This function reads the regrow file and store the values.
void PrepareRegrowFile(Parameters const & parameters, int ligandDescription[])
{
    // Open the regrow file, check if the file has been opened.
    ifstream regrowFile(parameters.regrowFile.c_str(), ios::in);
    if(!regrowFile) {
        cerr << "ERROR: The file " << parameters.regrowFile << " is missing, or I can't read it." << endl;
        exit (-1);   
        }

    // Store the information and close the file.
    for (int j=0 ; j < 3*parameters.maxFragment-2 ; j++) {
        regrowFile >> ligandDescription[j];
        }
    regrowFile.close();

    cout << endl << "*************** Reading of the regrow file successful ****************" << endl;
}

