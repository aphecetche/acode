
CXX := $(shell root-config --cxx)

CXXFLAGS := -g -O2 -Wall $(shell root-config --cflags)

CXXFLAGS += -I$(ALICE_ROOT)/include

LIBS := $(shell root-config --libs)

LIBS += -L$(ALICE_ROOT)/lib -L$(ALICE_PHYSICS)/lib

#LIBS += -lEG -lGeom -lVMC -lMinuit -lTree -lProof -lProofPlayer -lXMLParser -lPhysics -lSTEERBase -lESD -lAOD -lCDB -lRAWDatabase -lSTEER -lANALYSIS -lANALYSISalice -lHLTbase -lOADB

LIBS += -lSTEERBase -lESD -lAOD -lCDB -lRAWDatabase -lSTEER -lANALYSIS -lANALYSISalice -lHLTbase -lOADB -lProof -lPhysics -lEG

all: TestROOTFile rootFileSize rawDataTag dumpMC countEvents recPointMap branchSizes

countEvents: countEvents.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

rawDataTag: rawDataTag.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

rootFileSize: rootFileSize.o rootFileSizeMain.o
	$(CXX) $(CXXFLAGS) $^ $(shell root-config --libs) -o $@

branchSizes: branchSizes.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -L$(HOME)/o2/alfa/inst/lib -lboost_program_options -o $@

TestROOTFile: TestROOTFile.o TestROOTFileMain.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

dumpMC: dumpMC.o dumpMCmain.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

recPointMap: recPointMap.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -lMUONbase -lMUONrec -o $@

%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -I$(HOME)/o2/alfa/inst/include/boost -c $< -o $@

clean:
	rm -f *.o *.so *.d G__* TestROOTFile rootFileSize rawDataTag dumpMC countEvents recPointMap
