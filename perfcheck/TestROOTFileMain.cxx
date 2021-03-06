#include "TGrid.h"
#include "TestROOTFile.h"
#include <string>
#include <fstream>
#include "Riostream.h"
#include <map>
#include "TProof.h"
#include "TUrl.h"
#include "TFileInfo.h"
#include "TFileCollection.h"
#include "THashList.h"
#include "TObjectTable.h"

namespace {
const char* connect = "pod://";
}

Int_t TestDataSet(const char* dsname, const char* treename, std::map<std::string,int>& files)
{
  if (!gProof) TProof::Open(connect,"masteronly");
  if (!gProof) 
  {
    std::cout << "Could not connect to " << connect << std::endl;
    return -1;
    
  }
    
  std::cout << "Testing dataset " << dsname << " ... " << std::endl;
  
  TFileCollection* fc = gProof->GetDataSet(dsname);
  
  if (!fc)
  {
    std::cout << "cannot get dataset " << dsname << std::endl;
    return -2;
  }
  else
  {
    std::cout << "Scanning " << fc->GetList()->GetLast() << " files" << std::endl;
  }
  TIter next(fc->GetList());
  TFileInfo* fi;
  Int_t nbad(0);
  
  while ( ( fi = static_cast<TFileInfo*>(next()) ) )
  {
    TUrl url(*(fi->GetFirstUrl()));
    std::string filename(url.GetUrl());
    
//    TFileInfoMeta* meta = fi->GetMetaData();
//    
//    if (!meta)
//    {
////      continue; // no meta => not staged
//      std::cout << filename << " has no meta" << std::endl;
//    }
      if (!fi->TestBit(TFileInfo::kStaged) || fi->TestBit(TFileInfo::kCorrupted)) continue;
      
    int rv = TestROOTFile(filename.c_str(),treename);
    files[filename.c_str()] += rv;  
    if (rv<0) ++nbad;
  }
  
  std::cout << "nad=" << nbad << std::endl;
  
  return nbad;
}

int main(int argc, const char** argv)
{
  if ( argc < 2 )
  {
    std::cout << "usage " << argv[0] << " file(or dataset)name treename" << std::endl;
    return -1;
  }
  
  TString file(argv[1]);
    TString treename("aodTree");

  if ( argc > 2 )
  {
      treename = argv[2];
  }
  
  
  if (file.Contains(".root") && !file.BeginsWith("Find;") )
  {
    int rv = TestROOTFile(file.Data(),treename.Data());
    if (rv<0) std::cout << file.Data() << " : TestROOTFile FAILED" << std::endl;
    return 0;
  }
  else
  {
      std::vector<std::string> names;
      
    std::map<std::string,int> datasets;
    std::map<std::string,int> files;
    std::string dsname;

    Bool_t isDataSetList(kFALSE);
    
    if ( !file.Contains("#") && gSystem->AccessPathName(file.Data())==kFALSE) {
      
        std::ifstream in(file.Data());
      
        while ( in >> dsname )
        {
          // assume it's a text file containing a list of datasets to check
          
          if ( TString(dsname.c_str()).BeginsWith("Find;") ) {
            isDataSetList = kTRUE;
          }
          
          names.push_back(dsname);
        }
        in.close();
    }
    else
    {
        // assume it's the name of one single dataset
        names.push_back(file.Data());
    }
                        
    if (!isDataSetList)
    {
//      TString cmd;
//      cmd.Form("echo 0 $(awk '/Pss/ {print \"+\", $2}' /proc/%d/smaps) | bc",gSystem->GetPid());

      TStopwatch timer;
      
      for ( std::vector<std::string>::const_iterator it = names.begin(); it != names.end(); ++it )
      {
        int rv = TestROOTFile(it->c_str(),treename.Data());
        
//        gSystem->Exec(cmd.Data());

        if (rv<0) std::cout << it->c_str() << " : TestROOTFile FAILED" << std::endl;

//        gObjectTable->Print();
      }
      
      timer.Print();
      
      return 1;
    }

    // deal with the data set list
    
    for ( std::vector<std::string>::const_iterator it = names.begin(); it != names.end(); ++it )
    {
      std::string dsname = *it;
      datasets[dsname]=TestDataSet(dsname.c_str(),treename.Data(),files);
    }
    

    std::map<std::string,int>::const_iterator dsit;

    for ( dsit = datasets.begin(); dsit != datasets.end(); ++dsit ) 
    {
      Int_t n = dsit->second;
      if ( n ) 
      {
        std::cout << Form("Dataset %s has %d corrupted files",dsit->first.c_str(),dsit->second) << std::endl;
      }
    }
    
    std::map<std::string,int>::const_iterator it;
    for ( it = files.begin(); it != files.end(); ++it ) 
    {
      if ( it->second < 0 )
      {
        std::cout << Form("%3d %s",it->second,it->first.c_str()) << std::endl;
      }
    }
  }
  
  if ( gProof )
  {
    gProof->Close();
    gProof=0x0;
  }
  
  return 0;
}
