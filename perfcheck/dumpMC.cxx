#include "TGrid.h"
#include <string>
#include <fstream>
#include "Riostream.h"
#include <map>
#include "TTree.h"
#include "AliAODEvent.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "dumpMC.h"
#include "AliCDBManager.h"
#include "TSystem.h"
#include "AliRunLoader.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TKey.h"
#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"

//______________________________________________________________________________
dumpMC::dumpMC(int argc, char** argv) :
fFirstEvent(0),
fLastEvent(-1),
fOCDBPath("raw://"),
fPrimaryOnly(kTRUE),
fIsValid(kFALSE),
fFileName("")
{
  /// Main function for the program
  TObjArray args;
  
  for ( int i = 1; i < argc; ++i )
  {
    args.Add(new TObjString(argv[i]));
  }
  
  Int_t nok(0);
  
  for ( Int_t i = 0; i <= args.GetLast(); ++i )
  {
    TString a(static_cast<TObjString*>(args.At(i))->String());
    
    if ( !a.BeginsWith("--") && i == args.GetLast()  )
    {
      // last argument should be the filename
      continue;
    }
    
    if ( a == "--ocdb" )
    {
      fOCDBPath = static_cast<TObjString*>(args.At(i+1))->String();
      std::cout << "Using default storage  = " << fOCDBPath.Data() << std::endl;
      nok += 2;
      ++i;
    }
    else if( a == "--firstEvent")
    {
      fFirstEvent = static_cast<TObjString*>(args.At(i+1))->String().Atoll();
      nok+=2;
      ++i;
    }
    else if( a == "--lastEvent")
    {
      fLastEvent = static_cast<TObjString*>(args.At(i+1))->String().Atoll();
      nok+=2;
      ++i;
    }
    else if( a == "--all")
    {
      fPrimaryOnly = kFALSE;
      nok++;
    }
    else
    {
      if ( i != args.GetLast() )
      Usage();
    }
  }
  
  if ( nok < args.GetLast() - 1 )
  {
    Usage();
  }
  else
  {
    fFileName = static_cast<TObjString*>(args.At(args.GetLast()))->String();
  
//    AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
//    AliCDBManager::Instance()->SetRun(0);
    
    if (!fFileName.Contains(".root") )
    {
      Usage();
    }
    else
    {
      if ( fFileName.Contains("alien://") && !gGrid )
      {
        TGrid::Connect("alien://");
        if (!gGrid)
        {
          std::cerr << "cannot connect to the grid" << std::endl;
          fIsValid = kFALSE;
        }
      }
      fFileName = gSystem->ExpandPathName(fFileName.Data());
      fIsValid = kTRUE;
    }
  }
  
  std::cout << "filename=" << fFileName.Data() << std::endl;
}

//______________________________________________________________________________
Int_t dumpMC::Exec()
{
  if (!fIsValid) return -1;
  
  TFile* f = TFile::Open(fFileName.Data());
  if (f->IsZombie()) return -2;
  
  ULong64_t nevents(0);

  TTree* tree = static_cast<TTree*>(f->Get("aodTree"));
  if (tree)
  {
    nevents = tree->GetEntries();
    
    if (fLastEvent==-1)
    {
      fLastEvent = nevents-1;
    }
    
    fLastEvent = TMath::Min(fLastEvent,nevents);
    
    dumpAOD(tree);

    return 0;
  }

  // assume it's a kinematics
  
  TList* keys = f->GetListOfKeys();
  TIter next(keys);
  
  TKey* k;
  
  while ( ( k = static_cast<TKey*>(next()) ) )
  {
    TObject* object = k->ReadObj();
    
		if ( object->InheritsFrom("TDirectory") )
		{
      if (TString(object->GetName()).BeginsWith("Event")) ++nevents;
		}
  }
  
  if (!nevents)
  {
    std::cout << "Could not get events from file. " << fFileName.Data();
    std::cout << " Please check it is a Kinematics.root file" << std::endl;
    return -3;
  }
  
  if (fLastEvent==-1)
  {
    fLastEvent = nevents-1;
  }

  fLastEvent = TMath::Min(fLastEvent,nevents);
  
  dumpKine(*f);
  
  delete f;
  
  return 0;
}

//______________________________________________________________________________
void dumpMC::Print(const AliAODMCParticle& p) const
{
  TString name;
  
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(p.GetPdgCode());
  if (part)
  {
    name = part->GetName();
  }
  else
  {
    name = "Unknown";
    name += Form(" (pdgCode %d)",p.GetPdgCode());
  }

  std::cout << Form("%20s %10d Status %3d FirstMother %4d Daughters %4d %4d Eta %8.2f Pz %8.2f Pt %8.2f Vz %8.2f",
                    name.Data(),p.GetPdgCode(),p.GetStatus(),
                    p.GetMother(),
                    p.GetDaughter(0),
                    p.GetDaughter(1),p.Eta(),p.Pz(),p.Pt(),p.Zv()
                    ) << std::endl;
}

//______________________________________________________________________________
void dumpMC::Print(const TParticle& p) const
{
  std::cout << Form("%20s %10d Status %4d FirstMother %4d Daughters %4d %4d Eta %8.2f Pz %8.2f Pt %8.2f Vz %8.2f",
                    p.GetName(),p.GetPdgCode(),p.GetStatusCode(),
                    p.GetFirstMother(),p.GetFirstDaughter(),p.GetLastDaughter(),
                    p.Eta(),p.Pz(),p.Pt(),p.Vz()) << std::endl;
}

//______________________________________________________________________________
void dumpMC::GetMCGeneratorNames(const AliAODEvent& event)
{
  /// Find the generator(s) used for this MC
  
  TString geneNames;
  
  AliAODMCHeader* aodMCHeader = static_cast<AliAODMCHeader*>(event.FindListObject(AliAODMCHeader::StdBranchName()));
  
  if ( !aodMCHeader ) return;
  
  TList* lheaders = aodMCHeader->GetCocktailHeaders();
  AliGenEventHeader* gen;
  TIter next(lheaders);
  while ( ( gen = static_cast<AliGenEventHeader*>(next()) ) )
  {
    geneNames += gen->GetName();
    geneNames += " ";
  }

  std::cout << "---- File appear to be of generator(s) : " << std::endl;
  std::cout << geneNames.Data() << std::endl;
}

//______________________________________________________________________________
AliStack* dumpMC::GetStack(TTree* treeK)
{
  AliStack* stack = new AliStack(1000);
  
  TParticle* part(0x0);

  treeK->SetBranchAddress("Particles",&part);

//  Int_t n(0);
  
  for ( Int_t ip = 0; ip < treeK->GetEntries() ; ++ip )
  {
    treeK->GetEvent(ip);
    
    TVector3 pol;
    
    Print(*part);
    
    part->GetPolarisation(pol);
    
//    stack->PushTrack(0,
//                     part->GetFirstMother(),
//                     part->GetPdgCode(),
//                     part->Px(),
//                     part->Py(),
//                     part->Pz(),
//                     part->Energy(),
//                     part->Vx(),
//                     part->Vy(),
//                     part->Vz(),
//                     part->T(),
//                     pol.X(),
//                     pol.Y(),
//                     pol.Z(),
//                     kPNoProcess,
//                     n,
//                     1.0,
//                     part->GetStatusCode());
  }
  return stack;
}


//______________________________________________________________________________
void dumpMC::dumpKine(TFile& file)
{
//  Int_t entry;
//  if (id<fNprimary)
//    entry = id+fNtrack-fNprimary;
//  else
//    entry = id-fNprimary;
//  return entry;

  TParticle* part(0x0);
  
  for ( int i = fFirstEvent; i <= fLastEvent; ++i )
  {
    TTree* tree = static_cast<TTree*>(file.Get(Form("Event%d/TreeK",i)));

    if (!tree)
    {
      std::cout << "Cannot get Event" << i << "/TreeK" << std::endl;
      continue;
    }
    TObjArray particles(tree->GetEntries());
    particles.SetOwner(kTRUE);
    
    tree->SetBranchAddress("Particles",&part);
    
    std::cout << Form("----- Event %6d npart %6lld",i,tree->GetEntries()) << std::endl;

    Int_t firstPrimary(-1);
    
    for ( Int_t ip = 0; ip < tree->GetEntries() && firstPrimary < 0; ++ip )
    {
      tree->GetEntry(ip);
      if ( part->GetFirstMother()==-1 )
      {
        firstPrimary = ip;
      }
    }
    
//    std::cout << "firstPrimary is " << firstPrimary << std::endl;
    
    Int_t nprimaries = tree->GetEntries() - firstPrimary;
    
//    std::cout << "nprimaries is then " << nprimaries << std::endl;
    
    for ( Int_t ientry = 0; ientry < tree->GetEntries(); ++ientry )
    {
      tree->GetEntry(ientry);
      
      if (ientry < firstPrimary)
      {
        // secondaries go at the end
        particles.AddAt(new TParticle(*part),ientry+nprimaries);
      }
      else
      {
        particles.AddAt(new TParticle(*part),ientry-firstPrimary);
      }
    }
    
    TIter next(&particles);
    for ( Int_t ip = 0; ip < tree->GetEntries(); ++ip )
    {
      part = static_cast<TParticle*>(particles.At(ip));
      
      Bool_t mustShow(kTRUE);
      
      if ( fPrimaryOnly && ! (part->IsPrimary() && part->GetStatusCode()==1 ) )
      {
        mustShow = kFALSE;
      }
      
      if (mustShow)
      {
        std::cout << Form("%10d ",ip);
        if ( part )
        {
          Print(*part);
        }
        else
        {
          std::cout << " NULL " << std::endl;
        }
      }
    }
    
    delete part;
    part=0x0;
    
    std::cout << std::endl;
  }

}

//______________________________________________________________________________
void dumpMC::dumpAOD(TTree* tree)
{
  AliAODEvent event;
  event.ReadFromTree(tree);
  
  for ( ULong64_t i = fFirstEvent; i <= fLastEvent; ++i )
  {
    if ( tree->GetEntry(i) == 0 )
    {
      std::cout << "Read 0 bytes for entry " << i << std::endl;
      continue;      
    }

    TClonesArray *mcarray = static_cast<TClonesArray*>(event.FindListObject(AliAODMCParticle::StdBranchName()));

    if (!mcarray) continue;
    
    std::cout << Form("----- Event %6lld npart %6d",i,mcarray->GetLast()+1) << std::endl;
    
    TIter next(mcarray);
    AliAODMCParticle* mcPart;
    Int_t n(0);
    
    while ( ( mcPart = static_cast<AliAODMCParticle*>(next())) )
    {
      Bool_t mustShow(kTRUE);
      
      if ( fPrimaryOnly && ! (mcPart->IsPhysicalPrimary() && mcPart->GetStatus() ) )
      {
        mustShow = kFALSE;
      }
      
      if (mustShow)
      {
        std::cout << Form("%10d ",n);
        Print(*mcPart);
      }
      ++n;
    }
    
    std::cout << std::endl;
    
    if ( i == fLastEvent )
    {
      GetMCGeneratorNames(event);
    }
  }
}

//______________________________________________________________________________
void dumpMC::Usage()
{
  fIsValid = kFALSE;
  std::cout << "Usage : dumpMC [options] filename.root " << std::endl;
  std::cout << " where [options] is a combination of : " << std::endl;
  std::cout << "   --firstEvent n1 : first event to dump (default 0)" << std::endl;
  std::cout << "   --lastEvent n2 : last event to dump (default to last event on file) " << std::endl;
  std::cout << "   --ocdb ocdbPath : read the mapping from the given OCDB (default raw://)" << std::endl;
  std::cout << "   --all : show all particles (default is to show only primaries)" << std::endl;
}



