
CXX := $(shell root-config --cxx)

CXXFLAGS := -g -O2 -Wall $(shell root-config --cflags)

CXXFLAGS += -I$(ALICE_ROOT)/include

LIBS := $(shell root-config --libs)

LIBS += -L$(ALICE_ROOT)/lib

LIBS += -lGui -lEG -lGeom -lVMC -lMinuit -lTree -lProof -lProofPlayer -lXMLParser -lPhysics -lSTEERBase -lESD -lAOD -lCDB -lRAWDatabase -lSTEER -lANALYSIS -lANALYSISalice -lHLTbase -lSpectrum

all: rootstress

EventDict.cxx: Event.h EventLinkDef.h
	@echo "Generating dictionary $@..."
	rootcint -f $@ -c $^

rootstress: rootstress.o Event.o EventDict.o
	$(CXX) -Wall $(CXXFLAGS) $^ $(LIBS) -o $@

%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.so *.d G__* rootstress

