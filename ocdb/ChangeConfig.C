/*
 *  ChangeConfig
 *
 *  Get a config from OCDB and "patch" it
 *
 *
 */

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBRunRange.h"
#include "AliMUONCDB.h"
#include "AliMUON2DMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include <iostream>

AliMUON2DMap *PatchConfig(AliMUON2DMap &oldConfig,
                          const std::vector<int> &buspatches) {
  AliMUON2DMap *newConfig = new AliMUON2DMap(kTRUE);

  TIter next(oldConfig.CreateIterator());
  TObject* o;
  while ( (o = next())) {
      std::cout << o->ClassName() << "\n";
  }
  // for (auto bp : buspatches) {
  //   int deid = AliMpDDLStore::Instance()->GetDEFromBus(bp);
  //   if (!AliMpDEManager::IsValidDetElemId(deid)) {
  //     std::cerr << "Invalid de=" << deid << " from bp=" << bp << "\n";
  //     exit(1);
  //   }
  // }
  return newConfig;
}

//______________________________________________________________________________
void ChangeConfig(
    Int_t inputRun = 296935,
    const char *inputOCDB = "alien://folder=/alice/data/2018/OCDB",
    const char *outputOCDB =
        "alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2018") {
  AliCDBManager::Instance()->SetDefaultStorage(inputOCDB);
  AliCDBManager::Instance()->SetRun(inputRun);

  AliMpCDB::LoadAll();

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("MUON/Calib/Config");

  AliMUON2DMap *oldConfig = static_cast<AliMUON2DMap *>(entry->GetObject());

  AliCDBManager::Instance()->SetDefaultStorage(outputOCDB);

  auto bps = {1014, 1015, 1016, 1017};

  auto newConfig = PatchConfig(*oldConfig, bps);

  // AliMUONCDB::WriteToCDB(newConfig, "MUON/Calib/Config", 0,
  //                        AliCDBRunRange::Infinity(),
  //                        "Remove bus patches 1014, 1015, 1016, 1017 (see JIRA "
  //                        "ticket ALIROOT-8123)",
  //                        "L. Aphecetche");
}

