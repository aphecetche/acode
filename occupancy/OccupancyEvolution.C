#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliMUONTrackerConditionDataMaker.h"
#include "AliMUONVTrackerData.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "Riostream.h"
#include "TIterator.h"
#include "TMath.h"
#include "TObjString.h"
#include "TString.h"
#include <map>
#include <set>
#include <vector>
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONTrackerData.h"
#include "AliMergeableCollection.h"
#include "TH1F.h"
#include "AliMpConstants.h"
#include "TGraph.h"
#include "AliAnalysisMuMuGraphUtil.h"
#include "TCanvas.h"

///
/// Get the evolution of the occupancy, per station
///

//______________________________________________________________________________
void ReadIntegers(const char* filename, std::vector<int>& integers)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
  std::set<int> runset;
  
  char line[10000];
  
  in.getline(line,10000,'\n');
  
  TString sline(line);
  
  if (sline.Contains(","))
  {
    TObjArray* a = sline.Tokenize(",");
    TIter next(a);
    TObjString* s;
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }
    delete a;
  }
  else
  {
    runset.insert(sline.Atoi());
    
    while ( in >> i )
    {
      runset.insert(i);
    }
  }
  
  for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it ) 
  {
    integers.push_back((*it)); 
  }
  
  std::sort(integers.begin(),integers.end());
}

//______________________________________________________________________________
class RunInfo
{
public:
  RunInfo(Int_t run, 
          Long64_t cint1b, 
          Long64_t cmus1b,
          AliMUONVTrackerData* occupancy) : fNumber(run), fCINT1B(cint1b), fCMUS1B(cmus1b), 
  fOccupancy(occupancy) {
    
  }
  
  Int_t Number() const { return fNumber; }
  
  Long64_t NCINT1B() const { return fCINT1B; }

  Long64_t NCMUS1B() const { return fCMUS1B; }

  AliMUONVTrackerData* Occupancy() const { return fOccupancy; }

  Double_t EventSizePerDetectionElement(Int_t detElemId) const;

  Double_t EventSizePerChamber(Int_t chamberId) const;

  Double_t EventSizePerStation(Int_t stationId) const;

  void GetBusPatchesAboveOccupancy(std::set<int>& buspatches, Double_t occLimit=0.1) const;
  
private:
  Int_t fNumber;
  Long64_t fCINT1B;
  Long64_t fCMUS1B;
  AliMUONVTrackerData* fOccupancy;
};

//______________________________________________________________________________
Double_t RunInfo::EventSizePerDetectionElement(Int_t detElemId) const
{
  return fOccupancy->DetectionElement(detElemId,4); 
  // that's the number of hits in that detection elements for
  // fOccupancy->NumberOfEVents(ddl)
}

//______________________________________________________________________________
Double_t RunInfo::EventSizePerChamber(Int_t chamberId) const
{  
  AliMpDEIterator it;
  
  Double_t n(0);
  
  it.First(chamberId);
  
  while (!it.IsDone())
  {
    AliMpDetElement* de = it.CurrentDE();
    n += EventSizePerDetectionElement(de->GetId());
    it.Next();
  }
  return n;
}

//______________________________________________________________________________
Double_t RunInfo::EventSizePerStation(Int_t stationId) const
{
  return EventSizePerChamber(stationId*2) + EventSizePerChamber(stationId*2+1);
}

//______________________________________________________________________________
void RunInfo::GetBusPatchesAboveOccupancy(std::set<int>& buspatches, Double_t occLimit) const
{
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next()) ) )
  {
    Double_t occ = Occupancy()->BusPatch(bp->GetId(),2);
    
    if (occ>occLimit) buspatches.insert(bp->GetId());
  }
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& out, const RunInfo& run)
{
  Double_t meanOccupancy(0);
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;
  Int_t nbp(0);
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next()) ) )
  {
    Double_t occ = run.Occupancy()->BusPatch(bp->GetId(),2);

    meanOccupancy += occ;
    
    if (occ>0.1) cout << Form("BP %d occ %e",bp->GetId(),occ) << endl;
    
    ++nbp;
  }
  
  if (nbp) meanOccupancy /= nbp;
  
  cout << "BP <occ>=" << meanOccupancy << " nbp=" << nbp << endl;
  
  return out;
}

std::vector<RunInfo*> runs;

//______________________________________________________________________________
void OccupancyEvolutionBP(const char* runlist,
                          const char* ocdbPath="local:///Users/laurent/Alice/OCDBcopy2012")
{
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  
  Bool_t first(kTRUE);
  
  std::vector<int> runnumbers;
  
  ReadIntegers(runlist,runnumbers);

  for ( unsigned int i = 0 ; i < runnumbers.size(); ++i )
  {
    int runNumber = runnumbers[i];
    
    AliCDBManager::Instance()->SetRun(runNumber);

    if ( first ) 
    {
      AliMpCDB::LoadAll();  
//      AliMpCDB::LoadAll2();  
      first = kFALSE;
    }
    
    AliMUONTrackerConditionDataMaker occDM(runNumber, ocdbPath, "OCCUPANCY");
    
    occDM.SetOwnerOfData(kFALSE);
    
    AliMUONVTrackerData* occ = occDM.Data();
    
    runs.push_back(new RunInfo(runNumber,0,0,occ));    
  }
  
  std::vector<RunInfo*>::const_iterator it;
  
  for ( it = runs.begin(); it != runs.end(); ++it ) 
  {
    const RunInfo* ri = *it;
    if (!ri) continue;
    
    //if (ri) cout << (*ri) << endl;
    std::set<int> buspatches;
    ri->GetBusPatchesAboveOccupancy(buspatches,0.1);
    
    cout << Form("RUN %09d",ri->Number());
    
    for ( std::set<int>::const_iterator bit = buspatches.begin(); bit != buspatches.end(); ++bit )
    {
      cout << Form(" %4d",*bit);
    }
    cout << endl;
  }
}

//______________________________________________________________________________
AliMergeableCollection* GetManuOccupancyEvolution(const char* runlist,
                                               const char* ocdbPath="local:///Users/laurent/Alice/OCDBcopy2012")
{
  /// Get the manu occupancy histogrammed for all the run in runlist
  
  AliMergeableCollection* hc = new AliMergeableCollection("Occupancy","Occupancy histograms");
  
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  
  Bool_t first(kTRUE);
  
  std::vector<int> runnumbers;
  
  ReadIntegers(runlist,runnumbers);
  
//  AliCDBEntry* e = man->Get("MUON/Calib/OccupancyMap",runNumber);
//  
//  AliMUONVStore* occmap = static_cast<AliMUONVStore*>(e->GetObject());
//  
//  TIter next(occmap->CreateIterator());
//  AliMUONVCalibParam* p;
//  
//  while ( ( p = static_cast<AliMUONVCalibParam*>(next()) ) )
//  {
//    // convert buspatchid to detelemid
//    
//    Int_t busPatchId = p->ID0();
//    Int_t manuId = p->ID1();
//    
//    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);
//    
//    p->SetUniqueID(p->BuildUniqueID(detElemId,manuId));
//  }

  TH1* hB = new TH1F("BOccupancy","Global Bending Manu occupancy",1000,0,0.1);
  TH1* hNB = new TH1F("NBOccupancy","Global Non-Bending Manu occupancy",1000,0,0.1);

  hc->Adopt("/0000-AllRuns/",hB);
  hc->Adopt("/0000-AllRuns/",hNB);

  AliCDBManager* man = AliCDBManager::Instance();
  
  for ( unsigned int i = 0 ; i < runnumbers.size(); ++i )
  {
    int runNumber = runnumbers[i];

    std::cout << "Scanning run " << runNumber << std::endl;

    man->SetRun(runNumber);
    
    if ( first )
    {
      AliMpCDB::LoadAll();
      //      AliMpCDB::LoadAll2();
      first = kFALSE;
    }
    
    TH1* hoccB = new TH1F("manuBOccupancy","Bending Manu occupancy",1000,0,0.1);
    TH1* hoccNB = new TH1F("manuNBOccupancy","Non-Bending Manu occupancy",1000,0,0.1);
    
    hc->Adopt(Form("/%09d/",runNumber),hoccB);
    hc->Adopt(Form("/%09d/",runNumber),hoccNB);
    
    AliCDBEntry* e = man->Get("MUON/Calib/OccupancyMap",runNumber);
    AliCDBId id = e->GetId();
    if ( id.GetLastRun() == AliCDBRunRange::Infinity() )
    {
      std::cout << "ERROR : missing OccupancyMap for run " << runNumber << std::endl;
      continue;
    }
    
    AliMUONVStore* occmap = static_cast<AliMUONVStore*>(e->GetObject());
    TIter next(occmap->CreateIterator());
    AliMUONVCalibParam* p;
    while ( ( p = static_cast<AliMUONVCalibParam*>(next()) ) )
    {
//      Int_t detElemId = p->ID0();
      Int_t manuId = p->ID1();
      
      Double_t occ = p->ValueAsDouble(0,0) / p->ValueAsDouble(0,3) / p->ValueAsDouble(0,4);

      if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) )
      {
        hoccNB->Fill(occ);
        hNB->Fill(occ);
      }
      else
      {
        hoccB->Fill(occ);
        hB->Fill(occ);
      }
    }
  }
  
  return hc;
}

//______________________________________________________________________________
TGraph* GetFractionOfManuCutByOccupancyGraph(const AliMergeableCollection& hc, Double_t cut, Bool_t bending)
{
  /// Get, as a function of run number, the fraction of manu that have an occupancy above cut
  
  std::vector<Double_t> x;
  std::vector<Double_t> y;

  TList* list = hc.CreateListOfKeys(0);
  TIter next(list);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    Int_t runNumber;
    sscanf(str->String().Data(),"%09d",&runNumber);
    
    
    TH1* h = hc.Histo(str->String().Data(),Form("manu%sOccupancy",bending ? "B":"NB"));
    
    if ( !h ) continue;

    x.push_back(runNumber);

    Long_t n = 16828; //h->GetEntries();
    
    Int_t b = h->GetXaxis()->FindBin(cut);
    
    Double_t aboveCut = h->Integral(b,h->GetNbinsX());
    
    y.push_back(aboveCut/n);
  }

  
  delete list;
  
  return new TGraph(x.size(),&x[0],&y[0]);
}

//______________________________________________________________________________
Bool_t ManuOccupancyEvolution(const char* runlist,
                              const char* outputfile,
                              const char* ocdbPath="local:///Users/laurent/Alice/OCDBcopy2012")
{
  AliMergeableCollection* hc = GetManuOccupancyEvolution(runlist,ocdbPath);
  
  if (!hc)
  {
    std::cout << "ERROR : could not get manu occupancy for runlist" << std::endl;
    return kFALSE;
  }
  
  TFile* f = TFile::Open(outputfile,"RECREATE");
  
  if (!f || !f->IsOpen())
  {
    std::cout << "Cannot open outputfile " << outputfile << std::endl;
    delete hc;
    return kFALSE;
  }
  
  std::vector<double_t> cuts;
  
  cuts.push_back(0.01); // 1%
  cuts.push_back(0.015); // 1.5%
  cuts.push_back(0.03); // 3%
  
  AliAnalysisMuMuGraphUtil gu(ocdbPath);
  
  TObjArray graphsB;
  TObjArray graphsNB;
  
  for ( std::vector<Double_t>::size_type i = 0; i < cuts.size(); ++i )
  {
    std::cout << "Computing graph for cut " << cuts[i] << std::endl;
    
    for ( Int_t bending = 0; bending <= 1; ++bending )
    {
      TGraph* g = GetFractionOfManuCutByOccupancyGraph(*hc,cuts[i],bending);
      
      TString name;
      
      name.Form("FractionOf%sManuAbove%3.1fPercentOccupancy",
                bending ? "B":"NB",cuts[i]*100);
      
      name.ReplaceAll(".","_");
      
      if (!g) continue;
      

//      gu.Compact(*g);
      
      TGraph* gc = static_cast<TGraph*>(g->Clone(Form("occ(%s) > %3.1f %%",bending ? "Bending" : "Non Bending",cuts[i]*100)));
      
      if ( bending )
      {
        graphsB.Add(gc);
      }
      else
      {
        graphsNB.Add(gc);
      }

      g->Write(name.Data());
    }
  }

  TCanvas* cb = new TCanvas("bending","bending",800,400);
  TCanvas* cnb = new TCanvas("nonbending","nonbending",800,400);
  
  cb->cd();
  gu.PlotSameWithLegend(graphsB,0,0.04);
  
  cnb->cd();
  gu.PlotSameWithLegend(graphsNB,0,0.04);
  
  hc->Write();
  
  cb->Write();
  cnb->Write();
  
  f->Close();
  delete f;
  
  
  return kTRUE;
}

