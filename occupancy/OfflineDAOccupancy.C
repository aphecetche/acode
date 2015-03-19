///
/// Cut&Paste code to perform the same task as the MUONTRKOCCda, but reading from 
/// a .root file (instead of .raw file like the DA)

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

void OfflineDAOccupancy(Int_t runNumber, Int_t eventsPerRun=200);

void OfflineDAOccupancy(const char* file, Int_t eventsPerRun);

//______________________________________________________________________________
void Add(AliMUONVStore& destStore, const AliMUONVStore& srcStore)
{
  /// Add all elements from srcStore to destStore
  /// Each element of srcStore is supposed to be an AliMUONCalibParamNI,
  /// with ID0=busPatchId and ID1=manuId
  
  TIter next(srcStore.CreateIterator());
  AliMUONVCalibParam* source;
  
  while ( ( source = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    AliMUONCalibParamNI* dest = static_cast<AliMUONCalibParamNI*>(destStore.FindObject(source->ID0(),source->ID1()));
    if (!dest)
    {
      dest = static_cast<AliMUONCalibParamNI*>(source->Clone());
      destStore.Add(dest);
    }
    else
    {
      for ( Int_t i = 0; i < source->Size(); ++i ) 
      {
        for ( Int_t j = 0; j  < source->Dimension(); ++j ) 
        {
          dest->SetValueAsIntFast(i,j,dest->ValueAsIntFast(i,j)+source->ValueAsIntFast(i,j));
        }
      }
    }
  }
}

//______________________________________________________________________________
void GenerateOutputFile(const AliMUONVStore& store, ostream& out,
                        Int_t runNumber, Int_t nevents, Int_t numberOfEventsWithMCH)
{
  /// Write the channel hit count (grouped by manu) in the output file.
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* manu;
  
  out << "//===========================================================================" << endl;
  out << "//  Hit counter file calculated by " << __FILE__ << endl;
  out << "//===========================================================================" << endl;
  out << "//" << endl;
  out << "//       * Run Number          : " << runNumber << endl;
  out << "//       * File Creation Date  : " << TTimeStamp().AsString("l") << endl;
  out << "//---------------------------------------------------------------------------" << endl;
  out << "//  BP   MANU  SUM_N  NEVENTS" << endl;
  out << "//---------------------------------------------------------------------------" << endl;
  
  Int_t nlines(0);
  
  while ( ( manu = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t sum(0);
    
    for ( Int_t i = 0; i < manu->Size(); ++i )
    {
      sum += manu->ValueAsInt(i);
    }
    
    out << Form("%5d %5d %10d %10d",manu->ID0(),manu->ID1(),sum,nevents) << endl;
    ++nlines;
  }
  
  if (!nlines)
  {
    // ok : empty output. Is it because the run was empty ?
    if ( !numberOfEventsWithMCH )
    {
      // yes. Let give a hint in the file so the Shuttle preprocessor will know about this...
      out << Form("%5d %5d %10d %10d",-1,-1,0,0) << endl;
    }
    // else the preprocessor will fail, and that's good because there's something wrong somewhere...
  }
}

//______________________________________________________________________________
void OfflineDAOccupancy(Int_t runNumber, Int_t eventsPerRun)
{
  if (!gGrid)
  {
    TGrid::Connect("alien://");
  }
  if (!gGrid) return;
  
  const char* period="LHC10d";
  
  TGridCollection* gc = gGrid->OpenCollection(Form("alien:///alice/data/2010/%s/%09d/collection",period,runNumber));
  
  if ( !gc ) return;
  
  TFileCollection* fc = gc->GetFileCollection();
  
  TIter next(fc->GetList());
  TFileInfo* fi;
  TObjArray chunks;
  chunks.SetOwner(kTRUE);
  
  while ( ( fi = static_cast<TFileInfo*>(next()) ) )
  {
    TUrl* url = fi->GetFirstUrl();
    chunks.Add(new TObjString(Form("alien://%s",url->GetFile())));
  }

  int N = chunks.GetLast()/20;
  
  cout << "Will analyze " << N << " chunks" << endl;
  
  for ( int i = 0; i <= chunks.GetLast(); ++i )
  {
    TString file((static_cast<TObjString*>(chunks.At(i)))->String().Data());
    if ( ( i % N == 0 ) )
//    if (
//       file.Contains("23.130") 
//    || file.Contains("23.170") 
//    || file.Contains("23.240") ) 
    {
      OfflineDAOccupancy(file.Data(),eventsPerRun);
    }
  }
  
  delete fc;
  delete gc;
}

//______________________________________________________________________________
void OfflineDAOccupancy(const char* file, Int_t eventsPerRun)
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
  
  AliMUON2DMap oneEventData(kTRUE);
  AliMUON2DMap accumulatedData(kTRUE);

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
//    rawReader->SelectEquipment(-1,-1,-1);
//    rawReader->ReadHeader();
    
    oneEventData.Clear();
    
    ++numberOfEvents;
    
    if ( !rawReader->GetDataHeader() )
    {
      ++numberOfEventsWithoutCDH;
      
//      std::cout << Form("event %d type %d =? %d",numberOfEvents,rawReader->GetType(),
//                        rawReader->GetEventHeader()->Get("Type")) << std::endl;
    }
    
    if (rawReader->GetType() != AliRawEventHeaderBase::kPhysicsEvent)
    {
      if ( rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent )
      {
        ++numberOfCalibrationEvent;
      }
//      std::cout << Form("%p Skipping event of type %d",rawReader->GetDataHeader(),rawReader->GetType());
//      if ( rawReader->GetDataHeader() )
//      {
//        std::cout << Form("L1MSG=%x",rawReader->GetDataHeader()->GetL1TriggerMessage());
//      }
//      std::cout << std::endl;
//
      continue;
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
      AliMUONVCalibParam* one = static_cast<AliMUONVCalibParam*>(oneEventData.FindObject(buspatchId,manuId));
      
      if (!one)
      {
        one = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),buspatchId,manuId);
        oneEventData.Add(one);
      }
      
      one->SetValueAsInt(manuChannel,0,one->ValueAsInt(manuChannel,0)+1);
    }
    
    Bool_t badEvent = stream.HasPaddingError() || stream.HasGlitchError();
    
    if ( !badEvent )
    {
      ++numberOfUsedEvents;
      Add(accumulatedData,oneEventData);
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

  int year;
  int run;
  int chunk;
  int part;
  
  sscanf(gSystem->BaseName(file),"%02d%09d%03d.%d.root",&year,&run,&chunk,&part);
  
  cout << "Year " << year << " Run " << run << " Chunk " << chunk << " Part " << part << endl;
  
  ofstream fout(Form("occupancy.%09d%03d.%d",runNumber,chunk,part));
  
  GenerateOutputFile(accumulatedData,fout,runNumber,numberOfUsedEvents,numberOfEventsWithMCH);
  
  fout.close();
  
}
