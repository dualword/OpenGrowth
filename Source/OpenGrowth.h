#include <chrono>            // To initialize random numbers
#include <sys/types.h>       // To create directories
#include <sys/stat.h>        // To create directories
#include <random>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/atom.h>

const float PI = 3.14159265f;
const float BETA = 1.68780669f;        // BETA=1/RT. RT = (8.3144621*298.15)/(4.184*1000) = 0.59248491279 kcal/mol
const int MAX_FRAGMENTS = 500;

using namespace std;

struct Parameters{
    long double  seed;
    mt19937      random;
    string       rotamersName;
    string       conformersName;
    int          rotamersNumber;
    int          conformersNumber;
    float        proteinRange;
    string       mode;
    string       growthMode;
    string       probaFirstFrag;
    string       probaTransition;
    string       regrowFile;
    string       scoringFunction;
    string       energyFile;
    string       averageType;
    string       ligandName;
    double       bindingSite[3];
    float        bindingSize;
    string       fragmentListName;
    int          fragmentListSize;
    float        vdWScaleInter;
    float        vdWScaleIntra;
    float        branchingProba;
    int          rotationPrecision;
    int          optimizationMode;
    int          optimizationNumber;
    int          optimizationIterations;
    float        optimizationDistance;
    float        optimizationAngle;
    string       optimizationForceField;
    int          optimizationSteepestDescent;
    int          optimizationConjugateGradient;
    float        optimizationVdwCutOff;
    float        optimizationElecCutOff;
    float        MCTemp;
    int          maxFragment;
    unsigned int maxAtoms;
    int          maxMW;
    int          maxIterations;
    int          minFragments;
    unsigned int minAtoms;
    int          minEnergy;
    string       outputName;
    int          smilesOnly;
    int          writeDescription;
    string       threeMerFile;
    int          threeMerSize;
    int          verbose;
    int          numberOutputs;
    vector<unsigned int> forbiddenFragments;       // Number of the fragment in the array (start at 1).
};

struct Atom{
    double       coordinates[3];   // Coordinates of the atom.
    string       Type;             // Atom type for the scoring function when it is assigned as a string.
    int          TypeNumber;       // Atom type for the scoring function when it is assigned as a number.
    float        LJr;              // r0 value for the Lennard-Jones potential.
    float        LJe;              // Epsilon value for the Lennard-Jones potential.
    int          index;            // Atom number in the protein (we need it for StericClash).
};

struct Molecule{
    OpenBabel::OBMol     obMol;              // Molecule from OpenBabel.
    vector<unsigned int> growth;             // 0 for heavy atoms. For hydrogens, the type of fragment for this atom (i.e. place in the fragment array, start at 1).
    vector<int>          fragmentIndex;      // The index of the fragment in the molecule to which the atom belongs (6 if it was the 6th added fragment).
    vector<int>          fragmentNeighbours; // How many neighbours has the fragment to which this atom belongs.
    vector<unsigned int> growingSites;       // From which hydrogens can we grow. This array stores the index of the atom in the molecule, starting at 1.
    vector<Atom>         atom;               // To store the potential types of the atoms and their coordinates.
    int                  numberResidues;     // How many fragments in the molecule.
    float                interactionEnergy;  // Scoring energy.
    float                internalEnergy;     // Distortion energy (=energy(binding conformation)-energy(optimized geometry).
    int                  stericClashes;      // 0=none, 1=intra+inter, 2=only intra, 3=only inter.
    int                  ring;               // Used for fragments and not for ligands. 1=it is a ring.
    int                  aromatic;           // Used for correcting the protonation state. 1=it is an aromatic fragment.
};

struct Protein{
    OpenBabel::OBMol obMol;            // Molecule from OpenBabel.
    vector<Atom>     atom;             // To store the potential types of the protein atoms and their coordinates.
};

