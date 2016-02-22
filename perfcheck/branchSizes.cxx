
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <fstream>
#include <map>
#include "TTree.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TGrid.h"

namespace po = boost::program_options;
using namespace std;

typedef map<string,Long64_t> bs;


//_____________________________________________________________________________
void branchSizes(TTree& tree, bs& sizes, bool inMem)
{
  TObjArray* tb = tree.GetListOfBranches();

  if (!tb) return;

  tree.GetEntry(0);

  TBranch* b;
  TIter next(tb);

  while ( ( b = static_cast<TBranch*>(next()) ) )
  {
    if ( TString(b->GetName()).Contains("tracklets") ) continue;
        
    b->GetEntry(0);

    if ( inMem ) {
      sizes[b->GetName()] = b->GetTotBytes("*");
    } else {
      sizes[b->GetName()] = b->GetZipBytes("*");
    }
  }
}

//_____________________________________________________________________________
void rootFileSize(const char* filename, const char* treeName, bs& sizes, bool inMem)
{
  if (TString(filename).Contains("alien://"))
  {
    TGrid::Connect("alien://");
  }

  TFile* file = TFile::Open(filename);

  if (!file) return;

  TTree* tree = static_cast<TTree*>(file->Get(treeName));
  if (!tree) return;
  
  branchSizes(*tree,sizes,inMem);
  
  delete file;
}


//_____________________________________________________________________________
int main(int argc, char* argv[])
{
  string sFile;
  string sFileList;
  string sTreeName;
  vector<string> vFileList;
  bool inMem;
  
  try {
    
    po::options_description desc("Usage");
    desc.add_options()
    ("help,h", "produces this usage message")
    ("file,f", po::value<string>(&sFile), "name of a single file")
    ("filelist,l",po::value<string>(&sFileList), "name of a file list")
    ("tree,t",po::value<string>(&sTreeName),"name of the tree")
    ("mem,m",po::value<bool>(&inMem),"use the memory size and not the disk one")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (!vm.count("tree")) {
      sTreeName="aodTree";
    }
    
    if (vm.count("help")) {
      cout << desc << endl;
      return 1;
    }
    
    if (!sFile.empty()) {
      vFileList.push_back(sFile);
    }

    if (!sFileList.empty()) {
      string line;
      ifstream in(sFileList);
      while (getline(in,line))
      {
        vFileList.push_back(line);
      }
    }

  } catch (exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    cerr << "Exception of unknown type!\n";
  }

  bs sizes;
  map<string,bs> filebs;
  
  for ( vector<string>::size_type i =0 ; i < vFileList.size(); ++i ) {
    rootFileSize(vFileList[i].c_str(),sTreeName.c_str(),sizes,inMem);
    filebs[vFileList[i]]=sizes;
  }

  for ( map<string,bs>::const_iterator it = filebs.begin(); it != filebs.end(); ++it ) {
    string filename = it->first;
    cout << filename << endl;
    Long64_t total=0;
    for ( bs::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2 ) {
      total += it2->second;
    }
    for ( bs::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2 ) {
      cout << Form("%20c %20s %10lld (%5.1f%%)",' ',it2->first.c_str(),it2->second,100.0*it2->second/total) << endl;
    }
    cout << endl;
  }
    
  return 0;
}

