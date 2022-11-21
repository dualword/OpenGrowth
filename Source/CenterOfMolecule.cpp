#include <iostream>
#include <fstream>      // To write/read in files
#include <sstream>      // To convert int en string

using namespace std;

int main(int nbarg, char * argv[])
{
    cout << "*************************" << endl;
    cout << "**** Find the center ****" << endl;
    cout << "****  of a molecule  ****" << endl;
    cout << "*************************" << endl;
    cout << endl;

    ifstream molecule(argv[1], ios::in);
    if(!molecule) {
        cerr << "ERROR: The usage is ./CenterOfMolecule.exe File.xyz." << endl;
        exit (-1);
        }

    double Xmin=0, Xmax=0, Ymin=0, Ymax=0, Zmin=0, Zmax=0;

    string line;
    int moleculeSize;
    string moleculeName;
    int lineNumber=1;
    while(getline(molecule, line)) {
        istringstream input(line);
        string atom;
        float x, y, z;
        if      (lineNumber==1) { input >> moleculeSize; }
        else if (lineNumber==2) { input >> moleculeName; }
        else if (lineNumber==3) {
            input >> atom >> x >> y >> z;
            Xmin=x;
            Xmax=x;
            Ymin=y;
            Ymax=y;
            Zmin=z;
            Zmax=z;
            }
        else    {
            input >> atom >> x >> y >> z;
            if (x<Xmin) { Xmin=x; }
            if (x>Xmax) { Xmax=x; }
            if (y<Ymin) { Ymin=y; }
            if (y>Ymax) { Ymax=y; }
            if (z<Zmin) { Zmin=z; }
            if (z>Zmax) { Zmax=z; }
            }
        lineNumber++;
        }

    double Xcenter=0;
    double Ycenter=0;
    double Zcenter=0;
    Xcenter=(Xmin+Xmax)/2;
    Ycenter=(Ymin+Ymax)/2;
    Zcenter=(Zmin+Zmax)/2;

    cout << "Xmin    = " << Xmin    << "\tYmin    = " << Ymin    << "\tZmin    = " << Zmin    << endl;
    cout << "Xmax    = " << Xmax    << "\tYmax    = " << Ymax    << "\tZmax    = " << Zmax    << endl;
    cout << "Xcenter = " << Xcenter << "\tYcenter = " << Ycenter << "\tZcenter = " << Zcenter << endl;

    cout << endl;
    cout << "***************" << endl;
    cout << "****  End  ****" << endl;
    cout << "***************" << endl;
    cout << endl;

    return 0;
}

