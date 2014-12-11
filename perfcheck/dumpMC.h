#ifndef DUMPMC_H
#define DUMPMC_H

#include "TString.h"

class AliRunLoader;
class TParticle;
class AliAODMCParticle;
class TTree;
class TFile;
class AliStack;
class AliAODEvent;

class dumpMC
{
public:
  dumpMC(int argc, char** argv);
  
  Int_t Exec();
  
private:
  void Usage();
  void dumpAOD(TTree* tree);
  void dumpKine(TFile& file);
  void Print(const AliAODMCParticle& p) const;
  void Print(const TParticle& p) const;
  AliStack* GetStack(TTree* treeK);
  void GetMCGeneratorNames(const AliAODEvent& event);
  
private:
  ULong64_t fFirstEvent;
  ULong64_t fLastEvent;
  TString fOCDBPath;
  Bool_t fPrimaryOnly;
  Bool_t fIsValid;
  TString fFileName;
};

#endif
