#include "AliMUONVClusterStore.h"
#include "TFile.h"
#include "TKey.h"
#include "TList.h"
#include "Riostream.h"
#include "TMath.h"
#include "TTree.h"
#include "AliMergeableCollection.h"
#include "TH2F.h"
#include "AliMUONConstants.h"
#include <string>

void recPointMap(const char* file, AliMergeableCollection& hc)
{
  std::cout << file << "..." << std::endl;
  
  TFile* f = TFile::Open(file);
  
  if (!f) return;
  
  // get the number of events
  TList* keys = f->GetListOfKeys();
  TIter next(keys);
  
  TKey* k;
  int nevents(0);
  
  while ( ( k = static_cast<TKey*>(next()) ) )
  {
    TString event(k->GetName());
    event.ReplaceAll("Event","");
    nevents = TMath::Max(nevents,event.Atoi());
  }

  std::cout << "nevents=" << nevents << std::endl;

  AliMUONVClusterStore* clusterStore(0x0);
  
  for ( int i = 0; i < nevents; ++i )
  {
    TString object;
    object.Form("Event%d/TreeR",i);
    TTree* treeR = static_cast<TTree*>(f->Get(object.Data()));
    
    if (!treeR) continue;

    if ( !clusterStore ) 
    {
      clusterStore = AliMUONVClusterStore::Create(*treeR);
    }

    clusterStore->Clear();
    clusterStore->Connect(*treeR);
    treeR->GetEvent(0);

    TIter nextCluster(clusterStore->CreateIterator());
    AliMUONVCluster* cluster;
    
    while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster())) )
    {
      int ch = cluster->GetChamberId();
      TString hname;
      hname.Form("Chamber%d",ch);
      TH1* h = hc.Histo(hname.Data());

      if (!h)
      {
        Float_t rMax = AliMUONConstants::Rmax(ch/2);
        int nbins = 200; // 100
        
        h = new TH2F(hname, Form("cluster position distribution in chamber %d;X (cm);Y (cm)",ch+1), nbins, -rMax, rMax, nbins, -rMax, rMax);

        hc.Adopt(h);
      }

      h->Fill(cluster->GetX(),cluster->GetY());
    }
  }
  
  std::cout << "..." << nevents << " events treated" << std::endl;

  delete f;
}

int main(int argc, char** argv)
{
  if ( argc < 3 )
  {
    std::cout << "Usage : " << argv[0] << " inputfile(s) outputfile.root" << std::endl;
    return 1;
  }

  AliMergeableCollection hc("RP");

  TString inputName(argv[1]);

  TFile* f = TFile::Open(argv[2],"recreate");

  if (inputName.EndsWith(".root"))
  {
    recPointMap(argv[1],hc);
    hc.Write();
  }
  else
  {
    // assume it's a text file with the list of files to be treated
    ifstream in(inputName);
    std::string line;
    
    while (std::getline(in,line))
    {
      recPointMap(line.c_str(),hc);
      f->cd();
      hc.Write("",TObject::kOverwrite);
    }
    return 0;
  }
  
  delete f;
  
  return 0;
}