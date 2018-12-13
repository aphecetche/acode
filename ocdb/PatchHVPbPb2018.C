#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDCSValue.h"
#include "AliMUONCDB.h"
#include "AliMpDCSNamer.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

void KillChannel(const std::vector<std::string> &dcsNames, int runNumber,
                 const char *srcOCDBPath, const char *destOCDBPath);

void prepare(std::map<int, std::vector<std::string>> &run2hvtokill,
             const std::string &dcsName, const std::vector<int> &runs) {
  for (auto r : runs) {
    run2hvtokill[r].push_back(dcsName);
  }
}

void PatchHVPbPb2018() {
  const char *src = "alien://folder=/alice/data/2018/OCDB";
  const char *dest = "alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2018r";
  std::map<int, std::vector<std::string>> run2hvtokill;

  prepare(run2hvtokill, "MchHvLvRight/Chamber01Right/Quad4Sect3.actual.vMon",
          {297128});

  prepare(run2hvtokill, "MchHvLvLeft/Chamber03Left/Quad2Sect2.actual.vMon",
          {296694, 296749, 296836, 296890, 296899, 296966, 297193, 297218,
           297441, 297479, 297541, 297588});

  prepare(run2hvtokill, "MchHvLvLeft/Chamber04Left/Quad3Sect2.actual.vMon",
          {296694, 296750, 296849, 296890, 296899, 296930, 296966, 297315,
           297413, 297441, 297481, 297624});

  prepare(run2hvtokill, "MchHvLvLeft/Chamber04Left/Quad2Sect2.actual.vMon",
          {296794, 296935, 296979, 297222, 297446});

  prepare(run2hvtokill, "MchHvLvLeft/Chamber04Left/Quad3Sect1.actual.vMon",
          {296787, 297196, 297218, 297379});

  prepare(run2hvtokill, "MchHvLvLeft/Chamber03Left/Quad2Sect3.actual.vMon",
          {297031});

  prepare(run2hvtokill, "MchHvLvRight/Chamber04Right/Quad1Sect1.actual.vMon",
          {297372});

  prepare(run2hvtokill, "MchHvLvRight/Chamber04Right/Quad4Sect1.actual.vMon",
          {296932});

  int n{0};
  for (auto r : run2hvtokill) {
    n++;
    std::cout << std::setw(2) << n << " RUN " << r.first << " : killing ";
    for (auto s : r.second) {
      std::cout << s << ",";
    }
    std::cout << "\n";
    KillChannel(r.second, r.first, src, dest);
  }
}

void KillChannel(const std::vector<std::string> &dcsNames, int runNumber,
                 const char *srcOCDBPath, const char *destOCDBPath) {

  // function to patch the OCDB MUON/Calib/HV for HV channels that got into
  // trouble but were not correctly identified by the offline algorithm.
  // The dcs value for the HV channel (identified by its dcs _name_, not alias)
  // will be replaced by 0.

  AliCDBManager *man = AliCDBManager::Instance();

  AliMpDCSNamer dcsNamer("TRACKER");

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

    if (std::find(dcsNames.begin(), dcsNames.end(), name) != dcsNames.end()) {
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
  AliMUONCDB::WriteToCDB(hvmap, "MUON/Calib/HV", runNumber, runNumber,
                         "Patched for misbehaving channels (see ALIROOT-8155)",
                         "L. Aphecetche");
  man->ClearCache();
}
