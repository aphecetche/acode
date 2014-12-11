#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "TKey.h"
#include "Riostream.h"

//_____________________________________________________________________________
void branchSizes(TTree* tree, Long64_t& zipBytes, Long64_t& totBytes, TObjArray* lines)
{
  TObjArray* tb = tree->GetListOfBranches();

  if (!tb) return;

  tree->GetEntry(0);

  totBytes = tree->GetTotBytes();
  zipBytes = tree->GetZipBytes();

  lines->Add(new TObjString(Form("== Tree %20s TotalBytes %4.0f KB AfterCompression %4.0f KB Nentries %lld",
                                 tree->GetName(),totBytes/1024.0,zipBytes/1024.0,tree->GetEntries())));

  TBranch* b;
  TIter next(tb);
  
  Float_t check100(0.0);
  
  while ( ( b = static_cast<TBranch*>(next()) ) )
  { 
    b->GetEntry(0);

    TString msg(Form("      %20s TotalBytes %10d (%4.0f %%) AfterCompression %10d (%4.0f %%)",
                     b->GetName(),
                     (Int_t)b->GetTotBytes("*"),
                     (totBytes>0) ? b->GetTotBytes("*")*100.0/totBytes : 0,
                     (Int_t)b->GetZipBytes("*"),
                     (zipBytes>0) ? b->GetZipBytes("*")*100.0/zipBytes : 0
                     ));
    
    check100 += (zipBytes>0) ? b->GetZipBytes("*")*100.0/zipBytes : 0;
    
    lines->Add(new TObjString(msg));
  }

  TBranch* br = tree->BranchRef();

  if (br)
  {  
    TString msg(Form("      %20s AfterCompression %10d (%5.2f %% of file size)",
                     br->GetName(),(Int_t)br->GetZipBytes("*"),
                     br->GetZipBytes("*")*100.0/zipBytes
                     ));
    lines->Add(new TObjString(msg));
  }
  
  lines->Add(new TObjString(Form("check100 = %f",check100)));
  
}

//_____________________________________________________________________________
void Print(const TObjArray& lines)
{
  TIter next(&lines);
  TObjString* str;
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    if ( str->String().BeginsWith("==") ) std::cout << std::endl << std::endl;
    std::cout << str->String() << std::endl;
  }
}

//_____________________________________________________________________________
void rootFileSize(const char* filename, Bool_t showBranches)
{
  if (TString(filename).Contains("alien://"))
  {
    TGrid::Connect("alien://");
  }
  
  TFile* file = TFile::Open(filename);
  
  if (!file) return;

  Long64_t fileSize = file->GetSize();
  
  TIter nextKey(file->GetListOfKeys());
  TKey* key;
  
  TObjArray lines;
  lines.SetOwner(kTRUE);

  TObjArray branches;
  branches.SetOwner(kTRUE);

  while ( (key=static_cast<TKey*>(nextKey())) )
  {
    Long64_t diskSize = key->GetNbytes();
    Long64_t memSize = key->GetObjlen();
    
    if ( TString(key->GetClassName()) == "TTree" )
    {
      TTree* t = static_cast<TTree*>(key->ReadObj());
      branchSizes(t,diskSize,memSize,&branches);      
    }
    
    TString msg(Form("%20s size %10lld bytes on disk (in memory %10lld) (%5.2f %% of total file size)",
                      key->GetName(),
                      diskSize,
                      memSize,
                     100.0*diskSize/fileSize));
    
    lines.Add(new TObjString(msg));    
  }

  if ( showBranches ) 
  {
    Print(branches);
    std::cout << std::endl << std::endl;
  }
  
  Print(lines);
  
  delete file;
}
