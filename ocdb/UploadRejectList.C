#include "AliMpDCSNamer.h"
#include "AliMpDetElement.h"
#include "AliMpArrayI.h"
#include "AliCDBManager.h"
#include "AliMpCDB.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliMUONRejectList.h"
#include "Riostream.h"
#include "AliMpDDLStore.h"
#include "AliMpVSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpPCB.h"
#include "AliMpCathodType.h"
#include "AliMpSegmentation.h"
#include "AliMpMotifPosition.h"

void UploadRejectList(Int_t firstRun=0, Int_t lastRun=AliCDBRunRange::Infinity())
{
  AliCDBManager::Instance()->SetDefaultStorage("alien://folder=/alice/data/2015/OCDB");
  AliCDBManager::Instance()->SetRun(firstRun);
  AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/RejectList","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB2015");

  AliMpCDB::LoadAll();
  
  AliMUONRejectList rl;
  
  rl.SetDetectionElementProbability(806,1.0); // alignment problem
  
  // first round of checks in LHC13f with first run of the period : 4 holes
  rl.SetHVProbability("MchHvLvRight/Chamber03Right/Quad4Sect1",1.0);
  rl.SetPCBProbability(811,4,1.0);
  rl.SetPCBProbability(900,1,1.0);
  rl.SetPCBProbability(1010,0,1.0);
  
  // second scan of LHC13f: all runs
  
  // DE 612 : one of the 7-8 HV capa per PCB went bad ?
  // remove 3 full manus and then the channels in between...
  rl.SetManuProbability(612,1034,1.0);
  rl.SetManuProbability(612,112,1.0);
  rl.SetManuProbability(612,14,1.0);
  
  for ( Int_t i = 0; i < 8; ++i )
  {
    rl.SetChannelProbability(612,111,i,1.0);
    rl.SetChannelProbability(612,113,40+i,1.0);

    rl.SetChannelProbability(612,15,40+i,1.0);
    rl.SetChannelProbability(612,13,i,1.0);
    
  }
  
  for ( Int_t i = 48; i <= 63; ++i )
  {
    rl.SetChannelProbability(612,1141,i,1.0);
  }
  for ( Int_t i = 32; i <= 47; ++i )
  {
    rl.SetChannelProbability(612,1141,i,1.0);
  }
  
  for ( Int_t i = 0; i <= 31; ++i )
  {
    rl.SetChannelProbability(612,1140,i,1.0);
  }
  
  // DE 714 : same story as in DE 612 ...
  rl.SetManuProbability(714,111,1.0);
  rl.SetManuProbability(714,216,1.0);
  rl.SetManuProbability(714,1143,1.0);

  for ( Int_t i = 47; i <= 32; ++i )
  {
    rl.SetChannelProbability(714,112,i,1.0);
  }
  for ( Int_t i = 0; i <= 15; ++i )
  {
    rl.SetChannelProbability(714,215,i,1.0);
  }
  for ( Int_t i = 18; i <= 63; ++i )
  {
    rl.SetChannelProbability(714,1233,i,1.0);
  }
  for ( Int_t i = 32; i <= 47; ++i )
  {
    rl.SetChannelProbability(714,112,i,1.0);
  }
  rl.SetChannelProbability(714,1142,0,1.0);
  rl.SetChannelProbability(714,1142,1,1.0);

  
  rl.SetManuProbability(903,1129,1.0);
  rl.SetManuProbability(903,108,1.0);

  // only for run 196646
  rl.SetDetectionElementProbability(919,1.0);
  rl.SetPCBProbability(918,1,1.0);
  rl.SetPCBProbability(918,2,1.0);
  rl.SetPCBProbability(918,3,1.0);
  rl.SetBusPatchProbability(1648,1.0); // BP on DE 918
  
  AliCDBId id("MUON/Calib/RejectList",firstRun,lastRun);

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("MUON TRK");
  metaData.SetComment("Uploaded by UploadRejectList.C macro (handmade by L. Aphecetche) using information from MC vs Data tracking efficiency studies (Javier Martin Bianco) - eye scan (P. Pillot) - only for run 196646 -");
    
  AliCDBManager::Instance()->Put(&rl,id,&metaData);
}

