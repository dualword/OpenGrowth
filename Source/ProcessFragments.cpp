#include <iomanip>       // For setprecision
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/atom.h>

using namespace std;

// We define two structures which will be useful.
struct Atom{
    float x;
    float y;
    float z;
    string Type;
};

struct Molecule{
    Atom   atom[250];
    int    size;
    string name;
};

// A function to compute the distance between two atoms.
float dist(Atom const & atom1, Atom const & atom2)
{
    float Distance=sqrt( pow(atom1.x-atom2.x,2) + pow(atom1.y-atom2.y,2) + pow(atom1.z-atom2.z,2) );
    return Distance;
}

int main(int nbarg, char * argv[])
{
    // Write an error if not enough parameters.
    if (nbarg==1) {
        cerr << " " << endl;
        cerr << "ERROR: This program must be used either with \"./ProcessFragments.exe 1 File1.xyz\" or with" << endl;
        cerr << "\"./ProcessFragments.exe 2 SMARTS Library.smi\" where SMARTS is a SMARTS pattern. The first" << endl;
        cerr << "case will move the molecule File1.xyz, the second will count how many times the pattern is" << endl;
        cerr << "found in the Library.smi file." << endl;
        cerr << " " << endl;
        exit (-1);
    }

    // We will use the program for two things completely different. Either move a fragment or count SMARTS patterns in a library.
    string parameterEntry=argv[1];
    if (parameterEntry=="2")      {
        OpenBabel::OBConversion obConversion;
        OpenBabel::OBMol obMol;
        OpenBabel::OBFormat *format = obConversion.FormatFromExt(argv[3]);
        obConversion.SetInFormat(format);
        OpenBabel::OBSmartsPattern smarts;
        smarts.Init(argv[2]);
        vector<vector<int> > maplist;

        int total = 0;
        bool notatend = obConversion.ReadFile(&obMol, argv[3]);
        while (notatend) {
             smarts.Match(obMol);
             maplist = smarts.GetUMapList();
             total += maplist.size();
             obMol.Clear();
             notatend = obConversion.Read(&obMol);
             }
        cout << total << endl;
        }

    else if (parameterEntry=="1") {
    cout << "**************************************************************************" << endl;
    cout << "***  This program reads a .xyz file (optimized geometry), and puts the ***" << endl;
    cout << "***  first atom in the origin, the second in the x axis and the third  ***" << endl;
    cout << "***  in the xy plan. Use it with as many files as you want with:       ***" << endl;
    cout << "***  \"./ProcessFragments.exe 1 File1.xyz File2.xyz File3.xyz ...\"      ***" << endl;
    cout << "**************************************************************************" << endl;
    cout << endl;

    for(int i=2; i<nbarg; i++) {
        // Store the parameter names as strings.
        string rootName=argv[i];
        unsigned int positionToErase = rootName.find_last_of(".");
        rootName.erase(positionToErase);
        string fragmentFileName = rootName + ".xyz";
        string outputFileName   = rootName + "_New.xyz";
        cout << "Fragment to be processed: " << fragmentFileName << endl;

        // Open the fragment file, check if it has been opened.
        ifstream fragmentFile(fragmentFileName.c_str(), ios::in);
        if(!fragmentFile ) {
            cerr << "ERROR: The file " << fragmentFileName << " is missing, or I can't read it." << endl;
            exit (-1);
            }

        // Read the fragment file.
        Molecule fragment;
        string line;
        int lineCounter=1;
        while(getline(fragmentFile, line)) {
            istringstream input(line);
            if (lineCounter==1) {
                string sizeValue="";
                input >> sizeValue;
                fragment.size=atoi(sizeValue.c_str());
                }
            else if (lineCounter==2) {
                string nameValue="";
                input >> nameValue;
                fragment.name=nameValue;
                }
            else {
                string atomTypeValue="";
                string xValue="";
                string yValue="";
                string zValue="";
                input >> atomTypeValue >> xValue >> yValue >> zValue;
                fragment.atom[lineCounter-3].Type=atomTypeValue;
                fragment.atom[lineCounter-3].x=atof(xValue.c_str());
                fragment.atom[lineCounter-3].y=atof(yValue.c_str());
                fragment.atom[lineCounter-3].z=atof(zValue.c_str());
                }
            lineCounter++;
            }
        // Close the fragment file.
        fragmentFile.close();

        // Open the output file, and check if it has been opened.
        ofstream outputFile(outputFileName.c_str(), ios::out);
        if(!outputFile) {
            cerr << "ERROR: I cannot write in the output file " << outputFileName.c_str() << "." << endl;
            exit (-1);   
            }

        // Translation - The first atom of the first fragment is put in (0,0,0).
        Molecule translated;
        translated=fragment;
        for (int i=0 ; i < translated.size ; i++) {
            translated.atom[i].x = fragment.atom[i].x - fragment.atom[0].x;
            translated.atom[i].y = fragment.atom[i].y - fragment.atom[0].y;
            translated.atom[i].z = fragment.atom[i].z - fragment.atom[0].z;
            }

        // We define a fake atom for the x axis.
        Atom xAxisPoint;
        xAxisPoint.x = 1.0 ;
        xAxisPoint.y = 0.0 ;
        xAxisPoint.z = 0.0 ;

        // We will do a first rotation to align the 2 first atoms with the x axe. Let's define two vectors U (translated.atom[1]-translated.atom[0])
        // and V (xAxisPoint-translated.atom[0]) and normalize them. After the rotation, U will be aligned with V. We define W=UxV/||UxV||.
        // This operation is only made if U and V are not already aligned, i.e. if ||UxV||=0.
        Molecule rotated1;
        rotated1=translated;
        float uX = (translated.atom[1].x-translated.atom[0].x)/dist(translated.atom[1], translated.atom[0]);
        float uY = (translated.atom[1].y-translated.atom[0].y)/dist(translated.atom[1], translated.atom[0]);
        float uZ = (translated.atom[1].z-translated.atom[0].z)/dist(translated.atom[1], translated.atom[0]);
        float vX = (xAxisPoint.x-translated.atom[0].x)/dist(xAxisPoint, translated.atom[0]);
        float vY = (xAxisPoint.y-translated.atom[0].y)/dist(xAxisPoint, translated.atom[0]);
        float vZ = (xAxisPoint.z-translated.atom[0].z)/dist(xAxisPoint, translated.atom[0]);
        float normUV = sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
        float wX=0, wY=0, wZ=0;
        float cosAlpha=0, sinAlpha=0;
        if (normUV!=0) {
            wX = (uY*vZ - uZ*vY)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
            wY = (uZ*vX - uX*vZ)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));
            wZ = (uX*vY - uY*vX)/sqrt(pow(uY*vZ - uZ*vY,2) + pow(uZ*vX - uX*vZ,2) + pow(uX*vY - uY*vX,2));

            // We need to know from which angle we want to rotate. We know with the dot and the cross products: U.V = U*V*cos(alpha) and UxV = U*V*sin(alpha).
            // Since U and V are normalized to 1, the norm of UxV gives us |sin(alpha)|. This is not enough, we need to know the sign of alpha.
            // Thus, we use the following: sin(alpha) = det(W, U, V) where W=UxV/||UxV|| (the rule of Sarrus can be used). 
            cosAlpha = uX*vX + uY*vY + uZ*vZ;
            sinAlpha = (wX*uY*vZ + uX*vY*wZ + vX*wY*uZ) - (vX*uY*wZ + wX*vY*uZ + uX*wY*vZ);
            // Let's write the rotation matrix around the vector W (http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle).
            float rotationMatrix1[3][3];
            rotationMatrix1[0][0] = cosAlpha + pow(wX,2)*(1-cosAlpha);
            rotationMatrix1[0][1] = wX*wY*(1-cosAlpha) + wZ*sinAlpha;
            rotationMatrix1[0][2] = wX*wZ*(1-cosAlpha) - wY*sinAlpha;    // rotationMatrix1 has the following form:
            rotationMatrix1[1][0] = wX*wY*(1-cosAlpha) - wZ*sinAlpha;    // [0][0] [1][0] [2][0]
            rotationMatrix1[1][1] = cosAlpha + pow(wY,2)*(1-cosAlpha);   // [0][1] [1][1] [2][1]
            rotationMatrix1[1][2] = wY*wZ*(1-cosAlpha) + wX*sinAlpha;    // [0][2] [1][2] [2][2]
            rotationMatrix1[2][0] = wX*wZ*(1-cosAlpha) + wY*sinAlpha;
            rotationMatrix1[2][1] = wY*wZ*(1-cosAlpha) - wX*sinAlpha;
            rotationMatrix1[2][2] = cosAlpha + pow(wZ,2)*(1-cosAlpha);

            // Rotation 1 - The second atom is aligned with the xAxisPoint.
            for (int i=0 ; i < rotated1.size ; i++) {
                rotated1.atom[i].x = rotationMatrix1[0][0]*translated.atom[i].x + rotationMatrix1[1][0]*translated.atom[i].y + rotationMatrix1[2][0]*translated.atom[i].z ;
                rotated1.atom[i].y = rotationMatrix1[0][1]*translated.atom[i].x + rotationMatrix1[1][1]*translated.atom[i].y + rotationMatrix1[2][1]*translated.atom[i].z ;
                rotated1.atom[i].z = rotationMatrix1[0][2]*translated.atom[i].x + rotationMatrix1[1][2]*translated.atom[i].y + rotationMatrix1[2][2]*translated.atom[i].z ;
                }
            }

        // We will do a second rotation to put the third atom in the xy plan. We will rotate around the x axis from an angle Alpha.
        // This operation is only made if the third atom is not in the x axis.
        Molecule rotated2;
        rotated2=rotated1;
        wX = -1.0 ;
        wY = 0.0 ;
        wZ = 0.0 ;
        float distYZ=0;
        if(fragment.size>2) { distYZ = sqrt(pow(rotated1.atom[2].y,2) + pow(rotated1.atom[2].z,2)); }
        if(distYZ!=0) {
            cosAlpha = rotated1.atom[2].y / sqrt( pow(rotated1.atom[2].y,2) + pow(rotated1.atom[2].z,2) );
            sinAlpha = rotated1.atom[2].z / sqrt( pow(rotated1.atom[2].y,2) + pow(rotated1.atom[2].z,2) );
            float rotationMatrix2[3][3];
            rotationMatrix2[0][0] = cosAlpha + pow(wX,2)*(1-cosAlpha);
            rotationMatrix2[0][1] = wX*wY*(1-cosAlpha) + wZ*sinAlpha;
            rotationMatrix2[0][2] = wX*wZ*(1-cosAlpha) - wY*sinAlpha;    // rotationMatrix2 has the following form:
            rotationMatrix2[1][0] = wX*wY*(1-cosAlpha) - wZ*sinAlpha;    // [0][0] [1][0] [2][0]
            rotationMatrix2[1][1] = cosAlpha + pow(wY,2)*(1-cosAlpha);   // [0][1] [1][1] [2][1]
            rotationMatrix2[1][2] = wY*wZ*(1-cosAlpha) + wX*sinAlpha;    // [0][2] [1][2] [2][2]
            rotationMatrix2[2][0] = wX*wZ*(1-cosAlpha) + wY*sinAlpha;
            rotationMatrix2[2][1] = wY*wZ*(1-cosAlpha) - wX*sinAlpha;
            rotationMatrix2[2][2] = cosAlpha + pow(wZ,2)*(1-cosAlpha);

            // Rotation 2 - The third atom is put in the xy plan.
            for (int i=0 ; i < rotated2.size ; i++) {
                rotated2.atom[i].x = rotationMatrix2[0][0]*rotated1.atom[i].x + rotationMatrix2[1][0]*rotated1.atom[i].y + rotationMatrix2[2][0]*rotated1.atom[i].z ;
                rotated2.atom[i].y = rotationMatrix2[0][1]*rotated1.atom[i].x + rotationMatrix2[1][1]*rotated1.atom[i].y + rotationMatrix2[2][1]*rotated1.atom[i].z ;
                rotated2.atom[i].z = rotationMatrix2[0][2]*rotated1.atom[i].x + rotationMatrix2[1][2]*rotated1.atom[i].y + rotationMatrix2[2][2]*rotated1.atom[i].z ;
                }
            }

        // Write the output.
        outputFile << fixed << setprecision(6);
        if (fragment.size==2) {
            outputFile << rotated1.size << endl;
            outputFile << rotated1.name << endl;
            for (int i=0 ; i < rotated1.size ; i++) { outputFile << rotated1.atom[i].Type << "\t" << rotated1.atom[i].x << "\t" << rotated1.atom[i].y << "\t" << rotated1.atom[i].z << "\t0" << endl; }
            }
        else {
            outputFile << rotated2.size << endl;
            outputFile << rotated2.name << endl;
            for (int i=0 ; i < rotated2.size ; i++) { outputFile << rotated2.atom[i].Type << "\t" << rotated2.atom[i].x << "\t" << rotated2.atom[i].y << "\t" << rotated2.atom[i].z << "\t0" << endl; }
            }

        outputFile.close();
        cout << "\t" << outputFileName << " processed." << endl;
        }
    }
    return 0;
}

