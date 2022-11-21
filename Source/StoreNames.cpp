#include "OpenGrowth.h"

// This function stores the file names, since only one is given in the input file. From a name AAAA_0.pdb (which is read in the parameter value of ROTAMERS, CONFORMERS or LIGAND),
// it will write in the array given as a parameter: AAAA_0.pdb, AAAA_1.pdb, AAAA_2.pdb ... To do this, it removes everything after the last "_" character ("0.pdb"), and add "s.pdb"
// instead (where s is the number of the snapshot). For ligands, it works in the same way but with .xyz instead of .pdb.
void StoreNames(Parameters const & parameters, string snapshotNames[], string const typeOfSnapshots)
{
    // Define parameters.
    int snapshotNumber=0;
    string rootName="";
    string extension="";

    // Depending on which kind of snapshots we are dealing with, store the snapshotNumber and the rootName given in input.
    if (typeOfSnapshots == "rotamers") {
        cout << endl << "********** The following protein rotamers files will be used *********" << endl;
        snapshotNumber = parameters.rotamersNumber;
        rootName = parameters.rotamersName;
        extension = ".pdb";
        }
    else if (typeOfSnapshots == "conformers") {
        cout << endl << "********* The following protein conformers files will be used ********" << endl;
        snapshotNumber = parameters.conformersNumber;
        rootName = parameters.conformersName;
        extension = ".pdb";
        }
    else if (typeOfSnapshots == "ligandsRotamer") {
        cout << endl << "********** The following ligand files (rotamer) will be used *********" << endl;
        snapshotNumber = parameters.rotamersNumber;        // We assume the same number of ligands than rotamers.
        rootName = parameters.ligandName;
        extension = ".xyz";
        }
    else if (typeOfSnapshots == "ligandsConformer") {
        cout << endl << "********* The following ligand files (conformer) will be used ********" << endl;
        snapshotNumber = parameters.conformersNumber;      // We assume the same number of ligands than conformers.
        rootName = parameters.ligandName;
        extension = ".xyz";
        }
    // Keep only the beginning of the name.
    unsigned int positionToErase = rootName.find_last_of("_");
    rootName.erase(positionToErase);

    // Convert the name to string and then write the full name.
    for (int s=0 ; s < snapshotNumber ; s++) {
        // Counter.
        ostringstream oss_counter;
        oss_counter << s;
        string snapshotCounter = oss_counter.str();
        // Full snapshot name.
        string snapshotName = rootName + "_" + snapshotCounter + extension;
        snapshotNames[s]=snapshotName;
        cout << snapshotNames[s] << endl;
        }
    cout << "**********************************************************************" << endl;
}

