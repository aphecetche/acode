#include "TestROOTFile.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TFile.h"
#include "Riostream.h"
#include "TTree.h"

int ReadTree(const char* treename)
{
  TTree* tree = static_cast<TTree*>(gDirectory->Get(treename));
  if (!tree) return -2;
  
  Long64_t nentries = tree->GetEntries();
  
  for (Long64_t i = 0; i < nentries; ++i)
  {
    if ( tree->GetEntry(i) <= 0 ) return -2;
  }

  return nentries;
}

int TestROOTFile(const char* file, const char* treename)
{
  std::cout << "TestROOTFile " << file << " for tree " << treename << "..." << std::flush;
  Long64_t size(0);
  int rv(-1);
  
  if ( gSystem->AccessPathName(file) == 1 ) 
    {
      std::cout << " does not exists" << std::endl;
//      return 1;
    }
  else
    {
      TFile* f = TFile::Open(file);
      if (!f) 
      {
        std::cout << "Cannot open " << file << std::endl;
        return -1;
        
      }
        std::cout << " > " << std::flush;
      size += f->GetSize();
      if ( size > 0 ) 
      {
        rv = ReadTree(treename);
      }
        f->Close();
      delete f;
    }
  std::cout << Form("%10d entries read successfully",rv) << std::endl;
  
  if (rv<0) std::cout << "TestROOTFile : " << file << " has a problem" << std::endl;
  return rv;
}

