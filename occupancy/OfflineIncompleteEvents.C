
#include "AliMpConstants.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReaderDate.h"
#include "Riostream.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TGrid.h"
#include "TGridCollection.h"
#include "THashList.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TUrl.h"
#include "AliDAQ.h"

void OfflineIncompleteEvents(const char* file, Int_t eventsPerRun);

//______________________________________________________________________________
void OfflineIncompleteEvents(const char* file, Int_t eventsPerRun)
{
  
  if (!gGrid && TString(file).BeginsWith("alien://"))
  {
    TGrid::Connect("alien://");
    if (!gGrid) return;
  }
  
  AliRawReader* rawReader = AliRawReader::Create(file);
  
  int numberOfEvents(0);
  int numberOfPhysicsEvent(0);
  int numberOfCalibrationEvent(0);
  int numberOfEventsWithMCH(0);
  int numberOfEventsWithoutCDH(0);
  int numberOfIncompleteEvents(0);
  int numberOfFlushedEvents(0);

  int runNumber(-1);
  
  while (rawReader->NextEvent() && ( numberOfEvents < eventsPerRun || eventsPerRun < 0 ) )
  {
    rawReader->Reset();
    ++numberOfEvents;
    
    if ( !rawReader->GetDataHeader() )
    {
      ++numberOfEventsWithoutCDH;
    }
    
    if (rawReader->GetType() != AliRawEventHeaderBase::kPhysicsEvent)
    {
      if ( rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent )
      {
        ++numberOfCalibrationEvent;
        continue;
      }
    }
  
    if (runNumber<0) runNumber = rawReader->GetRunNumber();
    
    ++numberOfPhysicsEvent;
    
    if ( numberOfPhysicsEvent % 50 == 0 )
      cout << Form("%12d events processed : %12d physics %d calibration ones[ %d with MCH information ] and %d without CDH. %d incomplete. %d flushed.",
               numberOfEvents,numberOfPhysicsEvent,
               numberOfCalibrationEvent,
               numberOfEventsWithMCH,
               numberOfEventsWithoutCDH,
               numberOfIncompleteEvents,
               numberOfFlushedEvents) << endl;

    Bool_t mchThere(kFALSE);
    
    for ( int iDDL = 0; iDDL < AliDAQ::NumberOfDdls("MUONTRK") && !mchThere; ++iDDL )
    {
      rawReader->Reset();
      rawReader->Select("MUONTRK",iDDL,iDDL);
      if (rawReader->ReadHeader() )
      {
        if (rawReader->GetEquipmentSize() ) mchThere = kTRUE;
      }
    }
    
    if ( mchThere)
    {
      ++numberOfEventsWithMCH;
    }
    else
    {
      continue;
    }
  
  delete rawReader;
  
  cout << Form("%12d events processed : %12d physics %d calibration ones[ %d with MCH information ] and %d without CDH. %d incomplete. %d flushed.",
               numberOfEvents,numberOfPhysicsEvent,
               numberOfCalibrationEvent,
               numberOfEventsWithMCH,
               numberOfEventsWithoutCDH,
               numberOfIncompleteEvents,
               numberOfFlushedEvents) << endl;
}
