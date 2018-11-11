/*
 *  ChangeRecoParam.C
 *
 *  Created by Laurent Aphecetche on 23/08/10.
 *  
 *  Get the current valid recoparams from OCDB and "patch" them
 *
 * Last change 22-aug-2011 to lower CH3(R) HV limit
 *
 * Must be compiled and executed under aliroot (otherwise ITS libs not there)
 *
 */

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "TObjArray.h"
#include "AliDetectorRecoParam.h"
#include "Riostream.h"
#include "AliITSRecoParam.h"
#include "AliMUONRecoParam.h"
#include "AliMUONCDB.h"
#include "AliCDBRunRange.h"

//______________________________________________________________________________
void ChangeRecoParam(Int_t startRun=294307,const char* outputOCDB="local:///alice/data/2018/OCDB/")
{
  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetRun(99999999);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");
  
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
    if (AliRecoParam::Convert(p->GetEventSpecie())==AliRecoParam::kLowMult)
    {
      cout << Form("array[%d]=%s %s %s",i,
                   p ? p->ClassName() : "",
                   p ? AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(p->GetEventSpecie())) :"",
                   p ? ( p->IsDefault() ? "default" : "") : "" ) << endl;
      p->Print("");
      AliMUONRecoParam* rp = dynamic_cast<AliMUONRecoParam*>(p);
      if (!rp) 
      {
        cout << "OUPS. OUPS" << endl;
        return;
      }
      rp->SetPadGoodnessMask(0x8408BE9B);
      rp->Print("");
    }
  }
  
  AliCDBManager::Instance()->SetDefaultStorage(outputOCDB);
  
  AliMUONCDB::WriteToCDB(array, "MUON/Calib/RecoParam", startRun, AliCDBRunRange::Infinity(), 
                         "reconstruction parameters for MUON, patched for event size limits (for DQM)", "L. Aphecetche");
  
}

//______________________________________________________________________________
void ChangeITSRecoParam(Int_t startRun=157560,const char* outputOCDB="alien://folder=/alice/cern.ch/user/l/laphecet/OCDB")
{
  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetRun(startRun);
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("ITS/Calib/RecoParam");
  
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
    if (AliRecoParam::Convert(p->GetEventSpecie())==AliRecoParam::kLowMult)
    {
      cout << Form("array[%d]=%s %s %s",i,
                   p ? p->ClassName() : "",
                   p ? AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(p->GetEventSpecie())) :"",
                   p ? ( p->IsDefault() ? "default" : "") : "" ) << endl;
      AliITSRecoParam* rp = dynamic_cast<AliITSRecoParam*>(p);
      if (!rp) 
      {
        cout << "OUPS. OUPS" << endl;
        return;
      }
      rp->SetVertexer(4); // should be fast (i.e. MC) vertexer  (an enum would be better for sure)
    }
  }
  
  AliCDBManager::Instance()->SetDefaultStorage(outputOCDB);
  
  AliMUONCDB::WriteToCDB(array, "ITS/Calib/RecoParam", startRun, AliCDBRunRange::Infinity(), 
                         "reconstruction parameters for ITS, patched for fake vertex finder", "L. Aphecetche");
  
}


