#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

using namespace std;

//This program performs the 3MerScreen in a SMILES files. It works in the following way:
//    ./3MerScreen.exe File.smi Forbidden3Mers.dat

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

//////////////////////////////////////////////////////////////

int main(int nbarg, char * argv[])
{
    // Prepare the 3Mer list.
    int threeMerSize=SizeFile(argv[2]);
    string threeMerList[threeMerSize];
    ifstream threeMerFile(argv[2], ios::in);
    if(!threeMerFile) {
        cerr << "ERROR: The file " << argv[2] << " is missing, or I can't read it." << endl;
        exit (-1);   
        }
    // We read the file in argv[2] and store each line in the array threeMerList.
    string line;
    int j=0;
    while(getline(threeMerFile, line)) {
        threeMerList[j]=line;
        j++;
        }
    threeMerFile.close();

    // Read the SMILES file.
    OpenBabel::OBMol obMol;
    OpenBabel::OBConversion obConversion;
    obConversion.SetInFormat("smi");
    bool notatend = obConversion.ReadFile(&obMol, argv[1]);
    while (notatend) {
        int Counting=0;
        for (int z=0 ; z < threeMerSize && !Counting ; z++) {
            OpenBabel::OBSmartsPattern obSmarts;
            obSmarts.Init(threeMerList[z]);
            if (obSmarts.HasMatch(obMol)) { Counting++; }
            }
        if      (!Counting) { cout << obMol.GetTitle() << "\tPASS" << endl ; }
        else if (Counting)  { cout << obMol.GetTitle() << "\tFAIL" << endl ; }
        obMol.Clear();
        notatend = obConversion.Read(&obMol);
        }

    return 0;
}

