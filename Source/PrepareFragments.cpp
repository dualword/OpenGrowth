#include "OpenGrowth.h"

int SizeFile(string const fileName);

// This function prepares the fragments by storing their information (address, name, ring information, growth mode of atoms) in memory.
void PrepareFragments(Parameters & parameters, Molecule fragment[])
{
    parameters.fragmentListSize=SizeFile(parameters.fragmentListName);
    parameters.forbiddenFragments.clear();

    // To store the information from the fragment list.
    string address[parameters.fragmentListSize];

    // Open the fragment list file, check if the file has been opened.
    ifstream listFile(parameters.fragmentListName.c_str(), ios::in);
    if(!listFile) {
        cerr << "ERROR: The file " << parameters.fragmentListName << " is missing, or I can't read it." << endl;
        exit (-1);   
        }
    // Read each line of the fragment list file, store the information in address[] and close the file.
    string line;
    int j=0;
    while(getline(listFile, line)) {
        address[j]=line;
        j++;
        }
    listFile.close();

    // Store the rootName of the fragments (i.e. without "_1.xyz").
    string rootName[parameters.fragmentListSize];
    for (int i=0 ; i < parameters.fragmentListSize ; i++) {
        rootName[i]=address[i];
        unsigned int positionToErase = rootName[i].find_last_of("_");
        rootName[i].erase(positionToErase);
        }

    // For each fragment, store the information.
    for (int i=0 ; i < parameters.fragmentListSize ; i++) {
        // Create an object for OpenBabel.
        OpenBabel::OBConversion obConversion;
        OpenBabel::OBFormat *format = obConversion.FormatFromExt(address[i]);
        obConversion.SetInFormat(format);
        obConversion.ReadFile(&fragment[i].obMol, address[i]);

        // Open the fragment file, check if it has been opened.
        ifstream fragmentFile(address[i].c_str(), ios::in);
        if(!fragmentFile ) {
            cerr << "ERROR: The file " << address[i] << " is missing, or I can't read it." << endl;
            exit (-1);
            }

        // Prepare the fragment.
        fragment[i].growth.clear();
        fragment[i].fragmentIndex.clear();
        fragment[i].fragmentNeighbours.clear();
        int IsARing=0;
        int IsAromatic=0;

        // Check if this fragment is the same as one of the previous ones.
        int fragmentShift=0;
        for (int k=1 ; k <= i ; k++) {
            if (rootName[i]==rootName[i-k]) { fragmentShift++; }
            }

        // Read the fragment.
        unsigned int sizeFragment=0;
        string line;
        int lineCounter=1;
        while(getline(fragmentFile, line)) {
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
                // Update the dynamic arrays. If the growthMode is "random", we don't use the data read from the input file (atomGrowth) but store the index of the fragment in the fragment list instead.
                if (atomGrowth==0) { fragment[i].growth.push_back(0);                            }
                else               {
                    if (parameters.growthMode=="RANDOM") { fragment[i].growth.push_back(i+1);                        }
                    else                                 { fragment[i].growth.push_back(atomGrowth+i-fragmentShift); }
                    }
                fragment[i].fragmentIndex.push_back(1);
                fragment[i].fragmentNeighbours.push_back(0);
                }
            lineCounter++;
            }

        // Check if each atom is aromatic or is in a ring
        for (unsigned int j=1; j<=fragment[i].obMol.NumAtoms() ; j++ ) {
            // Get info from the atom.
            OpenBabel::OBAtom *atom;
            atom = fragment[i].obMol.GetAtom(j);
            if (fragment[i].growth[j-1]==1) {
                IsARing += atom->IsInRing();
                IsAromatic += atom->IsAromatic();
                }
            }

        // Define some additional parameters.
        fragment[i].growingSites.clear();
        for (unsigned int j=1; j<=fragment[i].growth.size() ; j++ ) {
            if(fragment[i].growth[j-1]>=1) { fragment[i].growingSites.push_back(j); }
            }
        fragment[i].numberResidues    = 1;
        fragment[i].interactionEnergy = 0;
        fragment[i].internalEnergy    = 0;
        fragment[i].stericClashes     = 0;
        if (IsARing == 0)    { fragment[i].ring = 0; }
        else                 { fragment[i].ring = 1; }
        if (IsAromatic == 0) { fragment[i].aromatic = 0; }
        else                 { fragment[i].aromatic = 1; }

        // Check that all molecule size are consistent, and close the file.
        if ( (fragment[i].obMol.NumAtoms()!=sizeFragment) || (fragment[i].obMol.NumAtoms()!=fragment[i].growth.size()) || (fragment[i].growth.size()!=sizeFragment) ) {
            cerr << "ERROR: some of the arrays describing the fragments have not the same size." << endl;
            exit (-1);
            }
        fragmentFile.close();
        }

    cout << endl << "************** Reading of fragment list file successful **************" << endl;
}

