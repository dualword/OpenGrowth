#include "OpenGrowth.h"

// This function prepares an array with the bad 3-mers by reading an input file.
void Prepare3Mer(Parameters const & parameters, string threeMerList[])
{
    if (parameters.threeMerFile != "") {
        // Open the file, check it can be read.
        ifstream threeMerFile(parameters.threeMerFile.c_str(), ios::in);
        if(!threeMerFile) {
            cerr << "ERROR: The file " << parameters.threeMerFile << " is missing, or I can't read it." << endl;
            exit (-1);   
            }

        // We read the file in parameters.threeMerFile and store each line in the array threeMerList.
        string line;
        int j=0;
        while(getline(threeMerFile, line)) {
            threeMerList[j]=line;
            j++;
            }
        threeMerFile.close();

        cout << endl << "**** Reading of the forbidden 3-mers successful (" << j << " patterns) ****" << endl;
        }
}

