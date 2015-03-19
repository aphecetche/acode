
#include "AliRawReader.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMergeableCollection.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "TDatime.h"
#include "AliMpCDB.h"
#include "TMath.h"
#include "AliRawEventHeaderBase.h"
#include "AliDAQ.h"
#include "TH1F.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"
#include "TProfile.h"
#include <vector>
#include <set>
#include "AliMpBusPatch.h"
#include "AliMpDDL.h"
#include "AliMpDEManager.h"

std::vector<int> timeResolutions;

//_________________________________________________________________________________________________
void GetTimeRange(Int_t runNumber, time_t& start, time_t& end)
{
  AliCDBManager::Instance()->SetRun(runNumber);

  AliCDBEntry* e = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  AliGRPObject* grp = static_cast<AliGRPObject*>(e->GetObject());

  start = grp->GetTimeStart();
  end = grp->GetTimeEnd();

  TDatime dstart(start);
  TDatime dend(end);

  dstart.Set(dstart.GetYear(),dstart.GetMonth(),dstart.GetDay(),dstart.GetHour(),0,0);
  dend.Set(dend.GetYear(),dend.GetMonth(),dend.GetDay(),dend.GetHour()+1,0,0);
  
  std::cout << TDatime(start).AsString() << std::endl;
  std::cout << TDatime(end).AsString() << std::endl;

  std::cout << dstart.AsString() << std::endl;
  std::cout << dend.AsString() << std::endl;

  start = dstart.Convert(kFALSE);
  end = dend.Convert(kFALSE);
}

//_________________________________________________________________________________________________
void FillCollection(AliMergeableCollection& hc, time_t runStart, time_t runEnd, int timeResolution)
{
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;
  
  Double_t xmin = runStart*1.0;
  Double_t xmax = runEnd*1.0;
  int nbins = TMath::Nint((runEnd-runStart)/timeResolution);
  
  // basis for all plot = per buspatch
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
  {
    TH1* h = new TH1F(Form("BP%04d",bp->GetId()),Form("Number of hits in %d s bins",timeResolution),nbins,xmin,xmax);
    
    h->GetXaxis()->SetTimeDisplay(1);
    h->GetXaxis()->SetTimeFormat("%H:%M");
    hc.Adopt(Form("/BUSPATCH/HITS/%ds",timeResolution),h);
  }
  
  nbins = TMath::Nint((runEnd-runStart)/timeResolution);
  
  // number of events needed for normalization
  
  TH1* h = new TH1F(Form("Nevents%ds",timeResolution),Form("Number of events %d s bins",timeResolution),nbins,xmin,xmax);
  
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetTimeFormat("%H:%M");
  
  hc.Adopt("",h);
}

//_________________________________________________________________________________________________
void FillNumberOfPads(AliMergeableCollection& hc)
{
  // number of pads needed for normalization

  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;

  Int_t total(0);
  
  AliMpDCSNamer dcs("TRACKER");
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
  {
    TH1* h = new TH1F(Form("BP%04d",bp->GetId()),"number of pads",1,0,1);
    
    Int_t npads(0);
    
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(bp->GetId());
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    for ( Int_t i = 0; i < bp->GetNofManus(); ++i )
    {
      Int_t manuId = bp->GetManuId(i);
      
      npads += de->NofChannelsInManu(manuId);
    }
    
    TH1* hde = hc.Histo(Form("/DE/NPADS/DE%04d",detElemId));
    if (!hde)
    {
        hde = new TH1F(Form("DE%04d",detElemId),"number of pads",1,0,1);
        hc.Adopt("/DE/NPADS/",hde);
    }
    hde->Fill(0.0,1.0*npads);
    
    Int_t ddlId = de->GetDdlId() + AliDAQ::DdlIDOffset("MUONTRK");
    
    TH1* hddl = hc.Histo(Form("/DDL/NPADS/DDL%d",ddlId));
    if (!hddl)
    {
      hddl = new TH1F(Form("DDL%d",ddlId),"number of pads",1,0,1);
      hc.Adopt("/DDL/NPADS/",hddl);
    }
    hddl->Fill(0.0,1.0*npads);
    
    Int_t chamberId = 1+AliMpDEManager::GetChamberId(detElemId);
    Int_t stationId = 1 + AliMpDEManager::GetChamberId(detElemId)/2;
    
    TH1* hchamberSide(0x0);
    TH1* hstationSide(0x0);
    
    if (dcs.DCSAliasName(detElemId).Contains("Left"))
    {
      hchamberSide = hc.Histo(Form("/CHAMBER/NPADS/CHAMBER%dLEFT",chamberId));
      if (!hchamberSide)
      {
        hchamberSide = new TH1F(Form("CHAMBER%dLEFT",chamberId),"number of pads",1,0,1);
        hc.Adopt("/CHAMBER/NPADS/",hchamberSide);
      }
      hstationSide = hc.Histo(Form("/STATION/NPADS/STATION%dLEFT",stationId));
      if (!hstationSide)
      {
        hstationSide = new TH1F(Form("STATION%dLEFT",stationId),"number of pads",1,0,1);
        hc.Adopt("/STATION/NPADS/",hstationSide);
      }
    }
    else
    {
      if (dcs.DCSAliasName(detElemId).Contains("Right"))
      {
        hchamberSide = hc.Histo(Form("/CHAMBER/NPADS/CHAMBER%dRIGHT",chamberId));
        if (!hchamberSide)
        {
          hchamberSide = new TH1F(Form("CHAMBER%dRIGHT",chamberId),"number of pads",1,0,1);
          hc.Adopt("/CHAMBER/NPADS/",hchamberSide);
        }
        hstationSide = hc.Histo(Form("/STATION/NPADS/STATION%dRIGHT",stationId));
        if (!hstationSide)
        {
          hstationSide = new TH1F(Form("STATION%dRIGHT",stationId),"number of pads",1,0,1);
          hc.Adopt("/STATION/NPADS/",hstationSide);
        }
      }
    }

    hchamberSide->Fill(0.0,1.0*npads);
    hstationSide->Fill(0.0,1.0*npads);
    
    TH1* hchamber = hc.Histo(Form("/CHAMBER/NPADS/CHAMBER%d",chamberId));
    if (!hchamber)
    {
      hchamber = new TH1F(Form("CHAMBER%d",chamberId),"number of pads",1,0,1);
      hc.Adopt("/CHAMBER/NPADS/",hchamber);
    }
    hchamber->Fill(0.0,1.0*npads);

    TH1* hstation = hc.Histo(Form("/STATION/NPADS/STATION%d",stationId));
    if (!hstation)
    {
      hstation = new TH1F(Form("STATION%d",stationId),"number of pads",1,0,1);
      hc.Adopt("/STATION/NPADS/",hstation);
    }
    hstation->Fill(0.0,1.0*npads);

    h->Fill(0.0,1.0*npads);
    
    total += npads;
    
    hc.Adopt("/BUSPATCH/NPADS/",h);
  }
  
  std::cout << "Total pads=" << total << std::endl;
  
}

//_________________________________________________________________________________________________
void GroupByStation(AliMergeableCollection& hc, int timeResolution)
{
  int station(1);
  
  for ( Int_t ich = 1; ich < 10; ich += 2 )
  {
    TH1* h = hc.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%d",timeResolution,ich));
    TH1* h1 = hc.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%d",timeResolution,ich+1));
    
    TH1* hstation = static_cast<TH1*>(h->Clone(Form("STATION%d",station)));
    
    hstation->Add(h1);
    hc.Adopt(Form("/STATION/HITS/%ds",timeResolution),hstation);

    h = hc.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dLEFT",timeResolution,ich));
    h1 = hc.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dLEFT",timeResolution,ich+1));
    
    hstation = static_cast<TH1*>(h->Clone(Form("STATION%dLEFT",station)));
    
    hstation->Add(h1);
    hc.Adopt(Form("/STATION/HITS/%ds",timeResolution),hstation);

    h = hc.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dRIGHT",timeResolution,ich));
    h1 = hc.Histo(Form("/CHAMBER/HITS/%ds/CHAMBER%dRIGHT",timeResolution,ich+1));
    
    hstation = static_cast<TH1*>(h->Clone(Form("STATION%dRIGHT",station)));
    
    hstation->Add(h1);
    hc.Adopt(Form("/STATION/HITS/%ds",timeResolution),hstation);
    
    ++station;
  }
}

//_________________________________________________________________________________________________
void GroupByChamber(AliMergeableCollection& hc, int timeResolution)
{
  for ( Int_t ich = 1; ich <= 10; ++ich )
  {
    AliMpDEIterator it;
  
    it.First(ich-1);
  
    TH1* hchamberLeft(0x0);
    TH1* hchamberRight(0x0);
    TList listLeft;
    TList listRight;
    listLeft.SetOwner(kFALSE);
    listRight.SetOwner(kFALSE);
    
    AliMpDCSNamer dcs("TRACKER");

    while (!it.IsDone())
    {
      Int_t detElemId = it.CurrentDEId();

       AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
      TH1* h = hc.Histo(Form("/DE/HITS/%ds/DE%04d",timeResolution,detElemId));
      
      if (dcs.DCSAliasName(detElemId).Contains("Left"))
      {
        if (!hchamberLeft)
        {
          hchamberLeft = static_cast<TH1*>(h->Clone(Form("CHAMBER%dLEFT",ich)));
        }
        else
        {
          listLeft.Add(h);
        }
      }
      else
      {
        if (!hchamberRight)
        {
          hchamberRight = static_cast<TH1*>(h->Clone(Form("CHAMBER%dRIGHT",ich)));
        }
        else
        {
          listRight.Add(h);
        }
      }
      
      it.Next();
    }
    
    hchamberLeft->Merge(&listLeft);
    hchamberRight->Merge(&listRight);
    hc.Adopt(Form("/CHAMBER/HITS/%ds",timeResolution),hchamberLeft);
    hc.Adopt(Form("/CHAMBER/HITS/%ds",timeResolution),hchamberRight);
    TH1* hchamber = static_cast<TH1*>(hchamberLeft->Clone(Form("CHAMBER%d",ich)));
    hchamber->Add(hchamberRight);
    hc.Adopt(Form("/CHAMBER/HITS/%ds",timeResolution),hchamber);
  }
}

//_________________________________________________________________________________________________
void GroupByDDL(AliMergeableCollection& hc, int timeResolution)
{
  Int_t nddls = AliDAQ::NumberOfDdls("MUONTRK");
  Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
  
  for ( Int_t i = 0; i < nddls; ++i )
  {
    Int_t ddlId = offset + i;
    
    AliMpDDL* ddl = AliMpDDLStore::Instance()->GetDDL(i);
    
    TH1* hddl(0x0);
    TList list;
    list.SetOwner(kFALSE);
    
    for ( Int_t ide = 0; ide < ddl->GetNofDEs(); ++ide )
    {
      Int_t detElemId = ddl->GetDEId(ide);
      
      TH1* h = hc.Histo(Form("/DE/HITS/%ds/DE%04d",timeResolution,detElemId));
      
      if (!hddl)
      {
        hddl = static_cast<TH1*>(h->Clone(Form("DDL%d",ddlId)));
      }
      else
      {
        list.Add(h);
      }
    }
    
    hddl->Merge(&list);
    hc.Adopt(Form("/DDL/HITS/%ds",timeResolution),hddl);
  }
}

//_________________________________________________________________________________________________
void GroupByDE(AliMergeableCollection& hc, int timeResolution)
{
  
  //  TH1* h = hc->Histo(Form("/PERBUSPATCH/HITS/%ds/%s",timeResolutions[is],bpName.Data()));

  //  AliMpDCSNamer dcs("TRACKER");
  //      TObjArray* a = dcs.DCSAliasName(detElemId).Tokenize("/");
  //
  //      TString chamberName = static_cast<TObjString*>(a->At(1))->String();
  //
  //      TString stationName(chamberName);
  //
  //      delete a;
  //
  //      int station(1);
  //
  //      for ( int i = 0; i < 10; i += 2 )
  //      {
  //        stationName.ReplaceAll(Form("Chamber%02d",i),Form("St%d",station));
  //        stationName.ReplaceAll(Form("Chamber%02d",i+1),Form("St%d",station));
  //        station++;
  //      }

  AliMpDEIterator it;
  
  it.First();
  
  while (!it.IsDone())
  {
    Int_t detElemId = it.CurrentDEId();
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    TList list;
    list.SetOwner(kFALSE);
    TH1* hde(0x0);
    
    if ( de->GetStationType() != AliMp::kStationTrigger)
    {
      for ( Int_t i = 0; i < de->GetNofBusPatches(); ++i )
      {
        Int_t busPatchId = de->GetBusPatchId(i);
        
        TH1* h = hc.Histo(Form("/BUSPATCH/HITS/%ds/BP%04d",timeResolution,busPatchId));
        
        if (!hde)
        {
          hde = static_cast<TH1*>(h->Clone());
          hde->SetName(Form("DE%04d",detElemId));
        }
        else
        {
          list.Add(h);
        }
      }
      
      hde->Merge(&list);
      hc.Adopt(Form("/DE/HITS/%ds",timeResolution),hde);
    }
    
    it.Next();
  }

}

//_________________________________________________________________________________________________
void OccupancyInTimeBins(const char* input, const char* output)
{
  timeResolutions.push_back(1);
  timeResolutions.push_back(10);
  timeResolutions.push_back(100);
  
  AliRawReader* rawReader = AliRawReader::Create(input);
  
  AliMUONRawStreamTrackerHP stream(rawReader);
  
  stream.DisableWarnings();
  stream.TryRecover(kTRUE);
  
  int numberOfUsedEvents(0);
  int numberOfBadEvents(0);
  int numberOfEvents(0);
  int numberOfPhysicsEvent(0);
  int numberOfCalibrationEvent(0);
  int numberOfEventsWithMCH(0);
  int numberOfEventsWithoutCDH(0);
  
  int runNumber(-1);
  
  time_t runStart, runEnd;
  AliMergeableCollection* hc(0x0);
  
  AliCDBManager* cdbm = AliCDBManager::Instance();
  
  if (!cdbm->IsDefaultStorageSet())
  {
    cdbm->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
  }
  
  cdbm->SetRun(0);
  
  AliMpCDB::LoadAll();


  while (rawReader->NextEvent() && numberOfEvents < 1000 )
  {
    rawReader->Reset();
    ++numberOfEvents;
    
    if ( !rawReader->GetDataHeader() )
    {
      ++numberOfEventsWithoutCDH;
    }
    
    if (rawReader->GetType() != AliRawEventHeaderBase::kPhysicsEvent)
    {
      if ( rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent )
      {
        ++numberOfCalibrationEvent;
      }
      continue;
    }
    
    if (runNumber<0)
    {
      runNumber = rawReader->GetRunNumber();
      GetTimeRange(runNumber,runStart,runEnd);

      hc = new AliMergeableCollection("occ");
      
      for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
      {
        FillCollection(*hc,runStart,runEnd,timeResolutions[is]);
      }
      
      FillNumberOfPads(*hc);
      
    }
    
    ++numberOfPhysicsEvent;
    
    if ( numberOfPhysicsEvent % 5000 == 0 )
      cout << Form("%12d events processed : %12d physics %d used ones %d bad ones [ %d with MCH information ]",
                   numberOfEvents,numberOfPhysicsEvent,numberOfUsedEvents,numberOfBadEvents,numberOfEventsWithMCH) << endl;
    
    Bool_t mchThere(kFALSE);
    
    for ( int iDDL = 0; iDDL < AliDAQ::NumberOfDdls("MUONTRK") && !mchThere; ++iDDL )
    {
      rawReader->Reset();
      rawReader->Select("MUONTRK",iDDL,iDDL);
      if (rawReader->ReadHeader() )
      {
        if (rawReader->GetEquipmentSize() ) mchThere = kTRUE;
      }
    }
    
    if ( mchThere)
    {
      ++numberOfEventsWithMCH;
    }
    else
    {
      continue;
    }
    
    Int_t buspatchId;
    UShort_t  manuId;
    UChar_t manuChannel;
    UShort_t adc;
    
    stream.First();
    
    std::map<int,int> bpValues;
    
    while ( stream.Next(buspatchId,manuId,manuChannel,adc,kTRUE) )
    {
      bpValues[buspatchId]++;
    }
    
    for ( std::map<int,int>::const_iterator it = bpValues.begin(); it != bpValues.end(); ++it )
    {
      const int& buspatchId = it->first;
      const int& bpvalue = it->second;
      
      TString bpName = Form("BP%04d",buspatchId);
      
      for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
      {
        TH1* h = hc->Histo(Form("/BUSPATCH/HITS/%ds/%s",timeResolutions[is],bpName.Data()));
        
        if (!h)
        {
          cout << "histogram not found" << endl;
          continue;
        }
        
        h->Fill(rawReader->GetTimestamp(),bpvalue);
      }
    }
    
    for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
    {
      TH1* h = hc->Histo(Form("Nevents%ds",timeResolutions[is]));
      
      if (!h)
      {
        cout << "histogram not found" << endl;
        continue;
      }
      
      h->Fill(rawReader->GetTimestamp());
    }
    
  }
  
  // Group BP histograms per DE then DDL then Chamber then Station
  
  for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
  {
    GroupByDE(*hc,timeResolutions[is]);
    GroupByDDL(*hc,timeResolutions[is]);
    GroupByChamber(*hc,timeResolutions[is]);
    GroupByStation(*hc,timeResolutions[is]);
  }
  
  // make normalized versions of the histograms
  TFile* fout = new TFile(output,"RECREATE");
  hc->Write("occ");
  delete fout;
}
