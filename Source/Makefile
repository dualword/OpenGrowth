
SRCGROW	=AddFragment.cpp CheckGrowth.cpp Energy.cpp Energy_SMOG2001.cpp Energy_SMOG2016.cpp FirstFragment.cpp FirstFragmentIfKnown.cpp LigandConstraints.cpp OpenGrowth.cpp OptimizationGeom.cpp OptimizationPosition.cpp Parse.cpp Prepare3Mer.cpp PrepareFragments.cpp PrepareProbabilityFiles.cpp PrepareProtein.cpp PrepareRegrowFile.cpp Random.cpp SaveOutput.cpp StericClash.cpp StoreNames.cpp
OBJGROW	=$(SRCGROW:.cpp=.o)

SRCFOG	=FOG2.cpp FOG2_AddFragment.cpp FOG2_Optimization.cpp FOG2_Parse.cpp FOG2_SaveOutput.cpp FOG2_StericClash.cpp
OBJFOG	=$(SRCFOG:.cpp=.o)

CC	=g++
CXXFLAGS =-O3 -Wall -std=c++11
LIBS    =-lm

### Debian, Ubuntu, Mint ###
INCLUDE =-I/usr/include/openbabel-2.0
LDFLAGS =-rdynamic /usr/lib/libopenbabel.so -Wl,-rpath,/usr/lib
### Fedora/OpenSuse ###
#INCLUDE =-I/usr/include/openbabel-2.0
#LDFLAGS =-rdynamic /usr/lib64/libopenbabel.so -Wl,-rpath,/usr/lib
### Local (OpenBabel compiled by yourself on Linux) ###
#INCLUDE =-I/home/nicolas/Programs/Openbabel_2.4.1/include/openbabel-2.0
#LDFLAGS =-rdynamic /home/nicolas/Programs/Openbabel_2.4.1/lib/libopenbabel.so -Wl,-rpath,/home/nicolas/Programs/Openbabel_2.4.1/lib
### MacOS (OpenBabel compiled by yourself on MacOS) ###
#INCLUDE =-I/Users/nicolas/Programs/Openbabel_2.4.1/include/openbabel-2.0
#LDFLAGS =/Users/nicolas/Programs/Openbabel_2.4.1/lib/libopenbabel.dylib -Wl,-rpath,/Users/nicolas/Programs/Openbabel_2.4.1/lib
### MinGW ###
#INCLUDE =-I"C:/Openbabel_2.4.1/include/openbabel-2.0"
#LDFLAGS =-L"C:/Openbabel_2.4.1/bin/" -lopenbabel.dll

all: OpenGrowth.exe FOG2.exe 3MerScreen.exe CenterOfMolecule.exe KBP2016-Training.exe ProcessFragments.exe SearchGUI.exe SMoG2016.exe clean

OpenGrowth.exe: $(OBJGROW)
	$(CC) $(OBJGROW) -o OpenGrowth.exe $(LDFLAGS)
FOG2.exe: $(OBJFOG)
	$(CC) $(OBJFOG) CheckGrowth.o Prepare3Mer.o PrepareFragments.o PrepareProbabilityFiles.o Random.o -o FOG2.exe $(LDFLAGS)
3MerScreen.exe: 3MerScreen.o
	$(CC) 3MerScreen.o -o 3MerScreen.exe $(LDFLAGS)
CenterOfMolecule.exe: CenterOfMolecule.o
	$(CC) CenterOfMolecule.o -o CenterOfMolecule.exe $(LDFLAGS)
KBP2016-Training.exe: KBP2016-Training.o
	$(CC) KBP2016-Training.o -o KBP2016-Training.exe $(LDFLAGS)
ProcessFragments.exe: ProcessFragments.o
	$(CC) ProcessFragments.o -o ProcessFragments.exe $(LDFLAGS)
SearchGUI.exe: SearchGUI.o
	$(CC) SearchGUI.o -o SearchGUI.exe $(LDFLAGS)
SMoG2016.exe: SMoG2016.o
	$(CC) SMoG2016.o -o SMoG2016.exe $(LDFLAGS)

.cpp.o:
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@
3MerScreen.o: 3MerScreen.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@
CenterOfMolecule.o: CenterOfMolecule.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@
KBP2016-Training.o: KBP2016-Training.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@
ProcessFragments.o: ProcessFragments.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@
SearchGUI.o: SearchGUI.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@
SMoG2016.o: SMoG2016.cpp
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@

clean:
	rm -f *.o

