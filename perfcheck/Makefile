
CXX := $(shell root-config --cxx)

CXXFLAGS := -g -Wall $(shell root-config --cflags)

#CXXFLAGS += -I$(ROOTSYS)/proof/proofd/inc

#CXXFLAGS += -I$(GSHELL_ROOT)/include/xrootd

CXXFLAGS += -I$(ALICE_ROOT)/RAW  -I$(ALICE_INSTALL)/include

LIBS := $(shell root-config --libs)

LIBS += -L$(ALICE_INSTALL)/lib/tgt_$(ALICE_TARGET)

LIBS += -lGui -lEG -lGeom -lVMC -lMinuit -lTree -lProof -lProofPlayer -lXMLParser -lPhysics -lSTEERBase -lESD -lAOD -lCDB -lRAWDatabase -lSTEER -lANALYSIS -lANALYSISalice -lHLTbase

#LIBS += -lXrdProofd -lProof -L$(GSHELL_ROOT)/lib -lXrdSys -lXrdNet

#LIBS += $(GSHELL_ROOT)/lib/libXrdOuc.a

all: TestROOTFile rootFileSize rawDataTag dumpMC countEvents

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

%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.so *.d G__* TestROOTFile rootFileSize rawDataTag dumpMC countEvents
