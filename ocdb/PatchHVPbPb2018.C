#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDCSValue.h"
#include "AliMUONCDB.h"
#include "AliMpDCSNamer.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <algorithm>
#include <iostream>
#include <vector>

void KillChannel(const std::vector<std::string> &dcsNames,
                 const std::vector<int> &runs, const char *srcOCDBPath,
                 const char *destOCDBPath);

void PatchHVPbPb2018() {
  const char *src = "alien://folder=/alice/data/2018/OCDB";
  const char *dest = "alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2018";

  KillChannel({"MchHvLvLeft/Chamber03Left/Quad2Sect2.actual.vMon"},
              {295854, 296068, 296195, 296241, 296269, 296511, 296548, 296619},
              src, dest);

  KillChannel({"MchHvLvLeft/Chamber04Left/Quad3Sect2.actual.vMon"},
              {295856, 295910}, src, dest);

  KillChannel({"MchHvLvLeft/Chamber04Left/Quad2Sect2.actual.vMon"},
              {296196, 296423}, src, dest);

  KillChannel({"MchHvLvLeft/Chamber04Left/Quad3Sect1.actual.vMon"},
              {296244, 296270, 296510, 296549}, src, dest);

  KillChannel({"MchHvLvLeft/Chamber04Left/Quad3Sect1.actual.vMon",
               "MchHvLvLeft/Chamber03Left/Quad2Sect2.actual.vMon"},
              {296304}, src, dest);
}

void KillChannel(const std::vector<std::string> &dcsNames,
                 const std::vector<int> &runs, const char *srcOCDBPath,
                 const char *destOCDBPath) {

  // function to patch the OCDB MUON/Calib/HV for HV channels that got into
  // trouble but were not correctly identified by the offline algorithm.
  // The dcs value for the HV channel (identified by its dcs _name_, not alias)
  // will be replaced by 0.

  AliCDBManager *man = AliCDBManager::Instance();

  AliMpDCSNamer dcsNamer("TRACKER");

  for (auto runNumber : runs) {

    man->SetDefaultStorage(srcOCDBPath);
    man->SetRun(runNumber);
    std::cout << "Run " << runNumber << std::endl;

    AliCDBEntry *entry = man->Get("MUON/Calib/HV");
    TMap *hvmap = static_cast<TMap *>(entry->GetObject()->Clone());

    TIter next(hvmap);
    TObjString *key;

    while ((key = static_cast<TObjString *>(next()))) {

      if (key->String().Contains("sw")) {
        continue;
      }

      auto name = dcsNamer.DCSNameFromAlias(key->String());

      if (std::find(dcsNames.begin(),dcsNames.end(),name) != dcsNames.end()) {
        std::cout << "Will kill channel name=" << name
                  << "\n(alias=" << key->String().Data() << ")"
                  << hvmap->GetValue(name) << "\n";
        TObjArray *a1 =
            static_cast<TObjArray *>(hvmap->GetValue(key->String().Data()));
        a1->Print();
        AliDCSValue *first = static_cast<AliDCSValue *>(a1->First());
        AliDCSValue *last = static_cast<AliDCSValue *>(a1->Last());
        auto start = first->GetTimeStamp();
        auto end = last->GetTimeStamp();
        a1->Delete();
        a1->Add(new AliDCSValue(0.0f, start));
        a1->Add(new AliDCSValue(0.0f, end));
        std::cout << " === \n";
        a1->Print();
      }
    }

    man->SetDefaultStorage(destOCDBPath);
    hvmap->SetUniqueID(hvmap->GetUniqueID() | (1 << 9));
    AliMUONCDB::WriteToCDB(
        hvmap, "MUON/Calib/HV", runNumber, runNumber,
        "Patched for misbehaving channels (see ALIROOT-8155)", "L. Aphecetche");
    man->ClearCache();
  }
}
