#include "AliRawReader.h"
#include "Riostream.h"
#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawVEquipment.h"
#include "AliRawEquipmentHeader.h"
#include "AliDAQ.h"

//    virtual AliRawEventHeaderBase *GetHeader() = 0;
//    virtual Int_t                  GetNEquipments() const = 0;
//    virtual AliRawVEquipment      *GetEquipment(Int_t index) const = 0;
//    virtual Int_t                  GetNSubEvents() const = 0;
//    virtual AliRawVEvent          *GetSubEvent(Int_t index) = 0;


void Test(const char* filename)
{
  AliRawReader* reader = AliRawReader::Create(filename);
  
  if (!reader) return;
  
  int nevent(0);
  
  while (reader->NextEvent())
  {
    cout << Form("event %5d detector pattern %p",nevent,reader->GetDetectorPattern()) << endl;
    ++nevent;
  }
}

void RawDataSanity(const char* filename)
{
  AliRawReader* reader = AliRawReader::Create(filename);
  
  if (!reader) return;
  
  int nevent(0);
  Double_t meanEventSize(0);
  
  while (reader->NextEvent())// && nevent < 1)
  {
    ++nevent;
    
    cout << Form("EVENT %d",nevent) << endl;
    
    AliRawVEvent* event = const_cast<AliRawVEvent*>(reader->GetEvent());
    
    Double_t eventSize = 0;
    
    for ( int i = 0; i < event->GetNSubEvents(); ++i ) 
    {
      AliRawVEvent* sub = event->GetSubEvent(i);

      for ( int j = 0; j < sub->GetNEquipments(); ++j ) 
      {
        AliRawVEquipment* eq = sub->GetEquipment(j);
        
        AliRawEquipmentHeader* equipmentHeader = eq->GetEquipmentHeader();

        UInt_t uid = equipmentHeader->GetId();
        
        int index;
        
        TString det(AliDAQ::DetectorNameFromDdlID(uid,index));
        
        if (det=="MUONTRK")     
        {
          cout << Form("%d %s %d",uid,det.Data(),equipmentHeader->GetEquipmentSize()) << endl;
          eventSize += equipmentHeader->GetEquipmentSize();
        }
      }
    }
    
    meanEventSize = meanEventSize*(1.0 - 1.0/nevent) + eventSize/nevent;
    
  }
  
  delete reader;
  
  cout << Form("Mean event size is %5.2f KB",meanEventSize/1024.0) << endl;
  
}