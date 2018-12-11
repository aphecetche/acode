/*
 *  ChangeConfig
 *
 *  Get a config from OCDB and "patch" it to remove some bus patches.
 *
 *  See https://alice.its.cern.ch/jira/browse/ALIROOT-8123
 *
 */

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBRunRange.h"
#include "AliMUON2DMap.h"
#include "AliMUONCDB.h"
#include "AliMUONCalibParamNF.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include <iostream>
#include <set>

// Count the number of buspatches present in the config.
int count(const AliMUON2DMap &config) {
  TIter next(config.CreateIterator());
  AliMUONVCalibParam *p;
  std::set<int> buspatches;
  while ((p = static_cast<AliMUONVCalibParam *>(next()))) {
    auto deid = p->ID0();
    auto manuid = p->ID1();
    buspatches.insert(AliMpDDLStore::Instance()->GetBusPatchId(deid, manuid));
  }
  return buspatches.size();
}

// Patch oldConfig by removing the buspatches present in the vector.
AliMUON2DMap *PatchConfig(AliMUON2DMap &oldConfig,
                          const std::vector<int> &buspatches) {
  AliMUON2DMap *newConfig = new AliMUON2DMap(kTRUE);

  TIter next(oldConfig.CreateIterator());
  AliMUONVCalibParam *p;
  while ((p = static_cast<AliMUONVCalibParam *>(next()))) {
    auto deid = p->ID0();
    auto manuid = p->ID1();
    auto bp = AliMpDDLStore::Instance()->GetBusPatchId(deid, manuid);
    if (std::find(buspatches.begin(), buspatches.end(), bp) !=
        buspatches.end()) {
      continue;
    }
    AliMUONVCalibParam *conf =
        static_cast<AliMUONVCalibParam *>(newConfig->FindObject(deid, manuid));
    if (!conf) {
      conf = new AliMUONCalibParamNF(1, 1, deid, manuid, 1);
      newConfig->Add(conf);
    }
  }
  return newConfig;
}

//______________________________________________________________________________
void ChangeConfig(Int_t inputRun, std::vector<int> bpToRemove,
                  const char *inputOCDB, const char *outputOCDB) {
  AliCDBManager::Instance()->SetDefaultStorage(inputOCDB);
  AliCDBManager::Instance()->SetRun(inputRun);

  AliMpCDB::LoadAll();

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("MUON/Calib/Config");

  AliMUON2DMap *oldConfig = static_cast<AliMUON2DMap *>(entry->GetObject());

  auto newConfig = PatchConfig(*oldConfig, bpToRemove);

  std::cout << "RUN " << inputRun << ": OLD CONFIG HAS " << count(*oldConfig)
            << " BUSPATCHES. NEW CONFIG HAS " << count(*newConfig)
            << " BUSPATCHES\n";

  AliCDBManager::Instance()->SetDefaultStorage(outputOCDB);
  AliMUONCDB::WriteToCDB(newConfig, "MUON/Calib/Config", inputRun, inputRun,
                         "Remove bus patches that were removed online but not "
                         "correctly propagated to OCDB "
                         "(see JIRA ticket ALIROOT-8123)",
                         "L. Aphecetche");
}

void ChangeConfigurations() {

  const char *inputOCDB = "alien://folder=/alice/data/2018/OCDB";
  const char *outputOCDB =
      "alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2018";

  for (auto run : {296935, 296938, 296941}) {
    ChangeConfig(run, {1014, 1015, 1016, 1017}, inputOCDB, outputOCDB);
  }

  for (auto run : {296967, 296968, 296969, 296971, 296972, 296973, 296974,
                   296975, 296976, 296977, 296979}) {
    ChangeConfig(run, {929,930}, inputOCDB, outputOCDB);
  }
}
