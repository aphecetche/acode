
CXX := $(shell root-config --cxx)

CXXFLAGS := -g -Wall $(shell root-config --cflags)

CXXFLAGS += -I$(ALICE_ROOT)/include

LIBS := $(shell root-config --libs)

LIBS += -L$(ALICE_ROOT)/lib

LIBS += -lGui -lEG -lGeom -lVMC -lMinuit -lTree -lProof -lProofPlayer -lXMLParser -lPhysics -lSTEERBase -lESD -lAOD -lCDB -lRAWDatabase -lSTEER -lANALYSIS -lANALYSISalice -lHLTbase

all: TestROOTFile rootFileSize rawDataTag dumpMC countEvents recPointMap

countEvents: countEvents.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

rawDataTag: rawDataTag.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

rootFileSize: rootFileSize.o rootFileSizeMain.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

TestROOTFile: TestROOTFile.o TestROOTFileMain.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

dumpMC: dumpMC.o dumpMCmain.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

recPointMap: recPointMap.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -lMUONbase -lMUONrec -o $@

%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.so *.d G__* TestROOTFile rootFileSize rawDataTag dumpMC countEvents recPointMap
