#include "TString.h"
#include "TObjArray.h"
#include "TTimeStamp.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "Riostream.h"
#include <vector>
#include <set>
#include <algorithm>
#include "CopyCDB.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMpDDLStore.h"
#include "AliMpBusPatch.h"
#include "AliMUONCDB.h"
#include "AliMpCDB.h"
#include "AliDetectorRecoParam.h"
#include "AliMUONRecoParam.h"

//______________________________________________________________________________
void ReadIntegers(const char* filename, std::vector<int>& integers)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  std::ifstream in(gSystem->ExpandPathName(filename));
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
void CopyCDB(const char* runlist, const char* fromURI, const char* toURI)
{
  std::vector<int> runs;
  
  ReadIntegers(runlist,runs);
  
  CopyCDB(runs,fromURI,toURI);
}

//______________________________________________________________________________
void CopyCDB(const std::vector<int>& runs, const char* fromURI, const char* toURI)
{
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i )
  {
    CopyCDB(runs[i],fromURI,toURI);
  }
}

//______________________________________________________________________________
void CopyCDB(Int_t runnr, const char* fromURI, const char* toURI)
{
  TString allowedObjects;

//  allowedObjects += " HLT/ConfigMUON/DecisionComponent";
//  allowedObjects += " HLT/ConfigMUON/FieldIntegrals";
//  allowedObjects += " HLT/ConfigMUON/HitReconstructor";
//  allowedObjects += " HLT/ConfigMUON/MansoTrackerFSM";
//  allowedObjects += " HLT/ConfigMUON/TriggerReconstructor";
  
  /*
   
  allowedObjects += " GRP/Geometry/Data";
  allowedObjects += " GRP/CTP/Config";
  allowedObjects += " GRP/GRP/LHCData";
  allowedObjects += " GRP/CTP/Scalers";
  allowedObjects += " GRP/GRP/Data";
  allowedObjects += " GRP/Calib/MeanVertexSPD";
  allowedObjects += " GRP/CTP/Aliases";

  allowedObjects += " MUON/Calib/GlobalTriggerCrateConfig";
  allowedObjects += " MUON/Calib/LocalTriggerBoardMasks";
  allowedObjects += " MUON/Calib/MappingData";
  allowedObjects += " MUON/Calib/RegionalTriggerConfig";
  allowedObjects += " MUON/Calib/TriggerLut";
   
  allowedObjects += " MUON/Calib/Config";
  allowedObjects += " MUON/Calib/Gains";
  allowedObjects += " MUON/Calib/GlobalTriggerCrateConfig";
  allowedObjects += " MUON/Calib/HV";
  allowedObjects += " MUON/Calib/LocalTriggerBoardMasks";
  allowedObjects += " MUON/Calib/MappingRunData";
  allowedObjects += " MUON/Calib/Neighbours";
  allowedObjects += " MUON/Calib/OccupancyMap";
  allowedObjects += " MUON/Calib/Pedestals";
  allowedObjects += " MUON/Calib/RecoParam";
  allowedObjects += " MUON/Calib/RegionalTriggerConfig";
  allowedObjects += " MUON/Calib/RejectList";
  allowedObjects += " MUON/Calib/TriggerDCS";
  allowedObjects += " MUON/Calib/TriggerEfficiency";
  allowedObjects += " MUON/Calib/TriggerLut";
  allowedObjects += " MUON/Calib/MappingData";

  allowedObjects += " MUON/Align/Data";

  allowedObjects += " GRP/Align/Data";
  allowedObjects += " ITS/Align/Data";
  allowedObjects += " VZERO/Align/Data";
  allowedObjects += " FMD/Align/Data";  
  allowedObjects += " T0/Align/Data";
  allowedObjects += " TPC/Align/Data";
  allowedObjects += " TRD/Align/Data";
  allowedObjects += " TOF/Align/Data";
  allowedObjects += " ACORDE/Align/Data";

  allowedObjects += " HLT/Calib/esdLayout";
  allowedObjects += " HLT/Calib/RecoParam";
  allowedObjects += " HLT/Calib/StreamerInfo";
  allowedObjects += " PHOS/Align/Data";
  allowedObjects += " EMCAL/Align/Data";
  allowedObjects += " HMPID/Align/Data";
   allowedObjects += " ZDC/Align/Data";
   allowedObjects += " PMD/Align/Data";
   
   
  
  allowedObjects += " GRP/Calib/MeanVertexTPC";
  allowedObjects += " GRP/Calib/CosmicTriggers";
  allowedObjects += " GRP/Calib/LHCClockPhase";
  allowedObjects += " GRP/CTP/CTPtiming";
  allowedObjects += " GRP/CTP/TimeAlign";
  allowedObjects += " GRP/Calib/RecoParam";
  allowedObjects += " GRP/CTP/Aliases";

  allowedObjects += " ITS/Calib/RecoParam";
  allowedObjects += " ITS/Calib/SPDNoisy";
  allowedObjects += " ITS/Calib/SPDDead";
  allowedObjects += " ITS/Calib/SPDSparseDead";
  allowedObjects += " ITS/Calib/CalibSDD";
  allowedObjects += " ITS/Calib/RespSDD";
  allowedObjects += " ITS/Calib/DriftSpeedSDD";
  allowedObjects += " ITS/Calib/DDLMapSDD";
  allowedObjects += " ITS/Calib/MapsTimeSDD";
  allowedObjects += " ITS/Calib/NoiseSSD";
  allowedObjects += " ITS/Calib/GainSSD";
  allowedObjects += " ITS/Calib/BadChannelsSSD";
  allowedObjects += " ITS/Calib/SPDFOEfficiency";
  
  allowedObjects += " ITS/Calib/SPDFONoise";
  
  allowedObjects += " TRIGGER/SPD/PITConditions";
  
  allowedObjects += " AD/Align/Data";
  */
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  // determine dynamically the current year
  TString fromUri(fromURI);
  cdb->SetDefaultStorage(fromUri.Data());
  cdb->SetRun(runnr);
  cdb->SetDrain(toURI);
  TString toUri(toURI);
  
  AliCDBStorage *defaultStorage = cdb->GetDefaultStorage();
  defaultStorage->QueryCDB(runnr);
  TObjArray* allIdsForRun = defaultStorage->GetQueryCDBList();
  TIter next(allIdsForRun);
  AliCDBId* id = 0;
  while ((id = dynamic_cast<AliCDBId*>(next())))
  {
    TString path(id->GetPath());
    if ( !allowedObjects.Contains(path.Data() ) ) continue;

    cdb->Get(path,cdb->GetRun());
  }
  
}

//______________________________________________________________________________
  void ChangeConfig(Int_t startRun, const char* fromURI,const char* toURI)
{
  AliCDBManager::Instance()->SetDefaultStorage(fromURI);
  AliCDBManager::Instance()->SetRun(startRun);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/Config");
  
  AliMpCDB::LoadAll();

  if (!entry) return;
  
  TObject* o = entry->GetObject();
  
  if (!o)
  {
    cout << "Could not get config ? Oups" << endl;
    return;
  }

  AliMUONVStore* config = static_cast<AliMUONVStore*>(o);

  AliMUONVStore* newConfig = static_cast<AliMUONVStore*>(config->Create());

  TIter next(config->CreateIterator());
  AliMUONVCalibParam* p;

  while ( ( p = static_cast<AliMUONVCalibParam*>(next())))
  {
    Int_t detElemId = p->ID0();
    Int_t manuId = p->ID1();

    Int_t bpid = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);

    AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(bpid);

    if ( bpid == 301 ) continue;
    if ( bpid == 302 ) continue;
    if ( bpid == 303 ) continue;
    if ( bpid == 313 ) continue;
    if ( bpid == 314 ) continue;
    if ( bpid == 315 ) continue;

    newConfig->Add(p->Clone());
  }

  AliCDBManager::Instance()->SetDefaultStorage(toURI);
  
  AliMUONCDB::WriteToCDB(newConfig, "MUON/Calib/Config", startRun, startRun,
                         "Config for MUON, patched for test removal of bus patches 301,302,303 L. Aphecetche");
  
}

//______________________________________________________________________________
void ChangeRecoParam(Int_t startRun,const char* fromURI,const char* toURI)
{
  AliCDBManager::Instance()->SetDefaultStorage(fromURI);
  AliCDBManager::Instance()->SetRun(startRun);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");
  
  AliMpCDB::LoadAll();

  if (!entry) return;
  
  TObject* o = entry->GetObject();
  
  if (!o)
  {
    cout << "Could not get recoparams ? Oups" << endl;
    return;
  }
  
  if ( o->IsA() != TObjArray::Class() ) 
  {
    cout << "This code only works with TObjArray recoparams. Sorry" << endl;
    return;
  }
  
  TObjArray* array = static_cast<TObjArray*>(o);
  for ( Int_t i = 0; i <= array->GetLast(); ++i ) 
  {
    AliDetectorRecoParam* p = static_cast<AliDetectorRecoParam*>(array->At(i));
    // if (AliRecoParam::Convert(p->GetEventSpecie())==AliRecoParam::kLowMult)
    // {
    //   cout << Form("array[%d]=%s %s %s",i,
    //                p ? p->ClassName() : "",
    //                p ? AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(p->GetEventSpecie())) :"",
    //                p ? ( p->IsDefault() ? "default" : "") : "" ) << endl;
    //   p->Print("");
      AliMUONRecoParam* rp = dynamic_cast<AliMUONRecoParam*>(p);
      if (!rp) 
      {
        cout << "OUPS. OUPS" << endl;
        return;
      }
//      rp->SetHVLimit(2,1580);
      UInt_t mask = rp->PadGoodnessMask();
      mask |= ( (1<<7) << 24 );
      rp->SetPadGoodnessMask(mask);
      rp->Print("");
  }
  
  AliCDBManager::Instance()->SetDefaultStorage(toURI);
  
  AliMUONCDB::WriteToCDB(array, "MUON/Calib/RecoParam", startRun, AliCDBRunRange::Infinity(), 
                         "reconstruction parameters for MUON, patched to take into account the bit for bus patch removed online by PAR", "L. Aphecetche");
  
}
