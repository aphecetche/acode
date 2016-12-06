#include "AliAnalysisTriggerScalers.h"
#include <iostream>
#include "AliCDBEntry.h"
#include "TMap.h"
#include "AliCDBManager.h"
#include "AliMUONCDB.h"
#include "AliDCSValue.h"
#include <vector>
#include <algorithm>

void PatchCDB(const char* runlist, const char* runlist1400, const char* srcOCDBPath="alien://folder=/alice/data/2016/OCDB", const char* destOCDBPath="alien://folder=/alice/cern.ch/user/l/laphecet/OCDBCH3L")
{
    // function to patch the OCDB MUON/Calib/HV for the swap of CH3L Q2S1 and Q2S2
    // runlist = full list of runs where Chamber03Left/Quad2Sect1 has an HV problem (trips, too low, plus the 1400 V
    // below)
    // runlist1400 = list of runs where Chamber03Left/Quad2Sect1 was struggling at 1400 V
    // for the runs in runlist1400, the HV will be forced to zero for that sector
    // note that Chamber03Left/Quad2Sect1 = Chamber02Left/Quad1Sect0 in DCS alias world
     
  AliAnalysisTriggerScalers ts(runlist,srcOCDBPath);

  std::vector<int> vrunlist = ts.GetRunList();

  AliAnalysisTriggerScalers ts1400(runlist1400,srcOCDBPath);
  std::vector<int> vrunlist1400 = ts1400.GetRunList();

  AliCDBManager* man = AliCDBManager::Instance();

  TObjString sector2("MchHvLvLeft/Chamber02Left/Quad1Sect0.actual.vMon");
  TObjString sector1("MchHvLvLeft/Chamber02Left/Quad1Sect1.actual.vMon");

  for ( auto r : vrunlist )
  {
      man->SetDefaultStorage(srcOCDBPath);
      man->SetRun(r);
      std::cout << "Run " << r << std::endl;

      AliCDBEntry* entry = man->Get("MUON/Calib/HV");
      TMap* hvmap = static_cast<TMap*>(entry->GetObject()->Clone());

      TPair* p1 = hvmap->RemoveEntry(&sector2);

      if ( std::find(vrunlist1400.begin(),vrunlist.end(),r) != vrunlist1400.end() )
      {
        TObjArray* a1 = static_cast<TObjArray*>(p1->Value());
        AliDCSValue* first = static_cast<AliDCSValue*>(a1->First());
        AliDCSValue* last = static_cast<AliDCSValue*>(a1->Last());
        a1->Delete();
        a1->Add(new AliDCSValue(0.0f,first->GetTimeStamp()));
        a1->Add(new AliDCSValue(0.0f,last->GetTimeStamp()));
      }
      TPair* p2 = hvmap->RemoveEntry(&sector1);

      hvmap->Add(new TObjString(sector2),p2->Value());
      hvmap->Add(new TObjString(sector1),p1->Value());

      delete p1->Key();
      delete p2->Key();

      man->SetDefaultStorage(destOCDBPath);
      hvmap->SetUniqueID( hvmap->GetUniqueID() | ( 1 << 9 ) );
      AliMUONCDB::WriteToCDB(hvmap,"MUON/Calib/HV",r,r,"Patched for CH3L Quad2Sect1 vs 0 swapping","L. Aphecetche");
      man->ClearCache();
  }
}
