
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONRawStreamTrackerHP.h"
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
#include "AliMUONLogger.h"

void OfflineReadoutErrors(const char* file, Int_t eventsPerRun);

//______________________________________________________________________________
void OfflineReadoutErrors(const char* file, Int_t eventsPerRun)
{
  
  if (!gGrid && TString(file).BeginsWith("alien://"))
  {
    TGrid::Connect("alien://");
    if (!gGrid) return;
  }
  
  AliRawReader* rawReader = AliRawReader::Create(file);
  
  AliMUONRawStreamTrackerHP stream(rawReader);
  
  stream.DisableWarnings();
  stream.TryRecover(kTRUE);

  AliMUONLogger* logger = new AliMUONLogger;

  if (logger)
  {
    stream.EnableMUONErrorLogger();  
    stream.SetMUONErrorLogger(logger);    
    stream.SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail);
  }

  int numberOfUsedEvents(0);
  int numberOfBadEvents(0);
  int numberOfEvents(0);
  int numberOfPhysicsEvent(0);
  int numberOfCalibrationEvent(0);
  int numberOfEventsWithMCH(0);
  int numberOfEventsWithoutCDH(0);
  
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
      cout << Form("%12d events processed : %12d physics %d used ones %d bad ones [ %d with MCH information ]",
                   numberOfEvents,numberOfPhysicsEvent,numberOfUsedEvents,numberOfBadEvents,numberOfEventsWithMCH) << endl;

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

    Int_t buspatchId;
    UShort_t  manuId;
    UChar_t manuChannel;
    UShort_t adc;
    
    stream.First();
    
    while ( stream.Next(buspatchId,manuId,manuChannel,adc,kTRUE) )
    {    
    }
    
    Bool_t badEvent = stream.HasPaddingError() || stream.HasGlitchError();
    
    if ( !badEvent )
    {
      ++numberOfUsedEvents;
    }
    else
    {
      ++numberOfBadEvents;
    }
    
  }
  
  std::cout << rawReader->ClassName() << std::endl;
  
  delete rawReader;
  
  cout << Form("%12d events processed : %12d physics %d used ones %d bad ones %d calibration ones[ %d with MCH information ] and %d without CDH",
               numberOfEvents,numberOfPhysicsEvent,numberOfUsedEvents,
               numberOfBadEvents,
               numberOfCalibrationEvent,
               numberOfEventsWithMCH,
               numberOfEventsWithoutCDH) << endl;

  logger->Print();
}
