///
/// Macro to change run range for a bunch of OCDB objects
/// in inputOCDB objects are supposed to be valid from 0 to infinity.

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "Riostream.h"

Bool_t ChangeRunRange(const char* objectPath,
                      int run1=0, int run2=AliCDBRunRange::Infinity(),
                      const char* inputOCDB="alien://folder=/alice/data/2013/OCDB",
                      const char* outputOCDB="alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2013") 
{
  AliCDBManager* man = AliCDBManager::Instance();
  
  man->SetDefaultStorage(inputOCDB);
  
  AliCDBEntry* e = man->Get(objectPath,AliCDBRunRange::Infinity());
  
  if (!e)
  {
    cout << Form("ERROR : could not get %s from %s",objectPath,inputOCDB) << endl;
    return kFALSE;
  }
  
  e->GetId().SetRunRange(run1,run2);
  
  AliCDBMetaData* md = e->GetMetaData();
  
  md->SetResponsible("L. Aphecetche and P. Pillot"); // to insure we have no $Id$ in the metadata fields (see https://savannah.cern.ch/bugs/?95527)
  
  man->SetDefaultStorage(outputOCDB);
 
  return man->Put(e->GetObject(),e->GetId(),e->GetMetaData());
}

