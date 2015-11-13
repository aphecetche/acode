/*
 *  ChangeRejectList
 *
 *  Get a rejectlist from OCDB and "patch" it
 *
 *
 */

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliMUONRejectList.h"
#include "AliMUONCDB.h"
#include "AliCDBRunRange.h"
#include "AliMpCDB.h"

//______________________________________________________________________________
void ChangeRejectList(Int_t inputRun=195861,
                      const char* inputOCDB="alien://folder=/alice/data/2013/OCDB",
                      const char* outputOCDB="alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2015")
{
  AliCDBManager::Instance()->SetDefaultStorage(inputOCDB);
  AliCDBManager::Instance()->SetRun(inputRun);
  
  AliMpCDB::LoadAll();

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RejectList");

  AliMUONRejectList* rl = static_cast<AliMUONRejectList*>(entry->GetObject());

  // repaired one
  rl->SetHVProbability("MchHvLvRight/Chamber03Right/Quad4Sect1",0.0);

  // new hole
  rl->SetPCBProbability(1023,3);

  AliCDBManager::Instance()->SetDefaultStorage(outputOCDB);
  
  AliMUONCDB::WriteToCDB(rl, "MUON/Calib/RejectList", 0, AliCDBRunRange::Infinity(),
                         "reject list for MUON, updated for 2015, from comparison of cluster maps between reconstructed data and simulations", "L. Aphecetche and P. Pillot");
}


