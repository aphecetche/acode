#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskRaw.h"
#include "AliRawInputHandler.h"
#include "AliRawReader.h"
#include "AliMUONTrackerDataMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrackerData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliLog.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"

ClassImp(AliAnalysisTaskRaw)

//________________________________________________________________________
AliAnalysisTaskRaw::AliAnalysisTaskRaw(const char *name, const char* ocdbpath) 
: AliAnalysisTask(name, ""), fDataMaker(0x0), fOCDBPath(ocdbpath)
{
  // Constructor

  AliInfo("");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a container
  DefineOutput(0, AliMUONTrackerData::Class());
}

//________________________________________________________________________
void AliAnalysisTaskRaw::ConnectInputData(Option_t *) 
{
  // Called once

  AliInfo("");

  AliRawInputHandler* eventHandler = dynamic_cast<AliRawInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  if (!eventHandler) 
  {
    Printf("ERROR: Could not get RawInputHandler");
  } 
}

//________________________________________________________________________
void AliAnalysisTaskRaw::CreateOutputObjects()
{
  AliInfo("");

}

//________________________________________________________________________
void AliAnalysisTaskRaw::LocalInit()
{
}

//________________________________________________________________________
void AliAnalysisTaskRaw::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliRawInputHandler* eventHandler = dynamic_cast<AliRawInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  if (!eventHandler) 
  {
    Printf("ERROR: Could not get RawInputHandler");
    return;
  } 
  
  AliRawReader* reader = eventHandler->GetRawReader();

  if (!AliMpDDLStore::Instance(kFALSE))
  {
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man->IsDefaultStorageSet())
    {
      man->SetDefaultStorage(fOCDBPath.Data());
      man->SetRun(reader->GetRunNumber());
    }
    else
    {
      AliWarning(Form("Default CDB storage was already set to %s. Could not set it to %s",
                      man->GetDefaultStorage()->GetURI().Data(),fOCDBPath.Data()));
    }
    AliMpCDB::LoadDDLStore();
    AliMpCDB::LoadManuStore();
  }
    
  if (!fDataMaker)
  {
      Bool_t histogram(kFALSE);
      AliMUONRecoParam* recoParam = AliMUONRecoParam::GetCosmicParam();
      
      AliInfo(Form("RUN=%d",reader->GetRunNumber()));
      
      fDataMaker = new AliMUONTrackerDataMaker(recoParam,
                                               reader->GetRunNumber(),
                                               reader,
                                               AliCDBManager::Instance()->GetDefaultStorage()->GetURI(),
                                               "NOGAIN",
                                               histogram,0.0,0.0);  
  }  
  
  fDataMaker->SetRawReader(reader);
  fDataMaker->ProcessEvent();
                             
  // Post output data.
  PostData(0, fDataMaker->Data());
}      

//________________________________________________________________________
void AliAnalysisTaskRaw::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  /*
  fHistPt = dynamic_cast<TH1F*> (GetOutputData(0));
  if (!fHistPt) {
    Printf("ERROR: fHistPt not available");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskRaw","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");
  */
}
