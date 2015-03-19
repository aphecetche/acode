#include "AliCDBManager.h"
#include "AliMpCDB.h"
#include "AliMUON2DMap.h"
#include "AliMUONTrackerIO.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include <fstream>
#include "Riostream.h"
#include "AliHLTDataTypes.h"

void ReadOccupancyFile(const char* filename, TString& out)
{
  std::ifstream in(filename);
  
  int headerSize = sizeof(AliHLTFXSHeader);
  
  in.ignore(headerSize);
  
  char line[1024];
  
  while ( in.getline(line,1023,'\n') )
  {
    out += line;
    out += "\n";
  }
}

void OccupancyAscii2OCDB(Int_t runNumber,
                         TObjArray& filenames,
                         const char* ocdbpath,
                         const char* comment)
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(runNumber);
  AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/OccupancyMap",ocdbpath);
  AliMpCDB::LoadAll();
  
  AliMUON2DMap occupancy(kTRUE);
  
	TObjString* str;
	TIter next(&filenames);
	while ( ( str = static_cast<TObjString*>(next()) ) )
	{
    TString data;
    ReadOccupancyFile(str->String().Data(),data);
    
	  AliMUONTrackerIO::DecodeOccupancy(data.Data(),occupancy);
	}
	
  if (occupancy.GetSize())
  {
    AliCDBId id("MUON/Calib/OccupancyMap",runNumber,runNumber);
    
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("MUON TRK");
    metaData.SetComment(comment);
    
    AliCDBManager::Instance()->Put(&occupancy,id,&metaData);
  }
}

void OccupancyAscii2OCDB(Int_t runNumber,
                         const char* filename,
                         const char* ocdbpath,
                         const char* comment)
{
  TObjArray a;
  a.SetOwner(kTRUE);
  a.Add(new TObjString(filename));
  OccupancyAscii2OCDB(runNumber,a,ocdbpath,comment);
}

void OccupancyAscii2OCDB(Int_t runNumber)
{
  OccupancyAscii2OCDB(runNumber,Form("occupancy.%d",runNumber),"alien://folder=/alice/cern.ch/user/l/laphecet/OCDB","Uploaded by Ascii2OCDB.C macro (using output from first occupancy DA)");
}

void OccupancyHLTtest()
{
//  tail -c +205 EOR_mch-occupancy_0x02_HLT\:FXS_CAL  | more > toto
// mv toto EOR_mch-occupancy_0x02_HLT\:FXS_CAL
  
  TObjArray hlt;
  hlt.SetOwner(kTRUE);
  
  hlt.Add(new TObjString("EOR_mch-occupancy_0x02_HLT:FXS_CAL"));
  hlt.Add(new TObjString("EOR_mch-occupancy_0x03_HLT:FXS_CAL"));
  
  OccupancyAscii2OCDB(196474,hlt,"alien://folder=/alice/cern.ch/user/l/laphecet/HCDB","HLT");
  OccupancyAscii2OCDB(196474,"occupancy.000196474082.15","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB",
                      "occupancy.000196474082.15");
}