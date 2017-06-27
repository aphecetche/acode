/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliRawInputHandler.h"

#include "AliLog.h"
#include "AliRawVEvent.h"
#include "AliRawReaderRoot.h"
#include "TTree.h"
#include "AliRawEventHeaderBase.h"

ClassImp(AliRawInputHandler)

//_____________________________________________________________________________
AliRawInputHandler::AliRawInputHandler() : 
AliVEventHandler(), 
fTree(0x0), 
fRawReader(0x0),
fRawEvent(0x0)
{
}

//_____________________________________________________________________________
AliRawInputHandler::AliRawInputHandler(const char* name, const char* title) : 
AliVEventHandler(name,title), 
fTree(0x0), 
fRawReader(0x0),
fRawEvent(0x0)
{
}

//_____________________________________________________________________________
AliRawInputHandler::~AliRawInputHandler()
{
  delete fRawEvent;
  delete fRawReader;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::Init(Option_t* /*opt*/)
{
  /// opt can be "proof" or "local"
//  AliInfo(Form("opt=%s",opt));
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::Init(TTree* tree, Option_t*)
{
//  AliInfo("");
  if (!tree)
  {
    AliError("got a null tree");
    return kFALSE;
  }
  fTree = tree;
//  fTree->Print();
//  AliInfo(Form("fTree=%p",fTree));
  
  delete fRawEvent;
  delete fRawReader;
  
  fRawEvent = 0x0;
  fRawReader = new AliRawReaderRoot; // must be created here, as this
  // object must be "fully" fonctional after Init, e.g. the GetRawReader
  // should return something...
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::BeginEvent(Long64_t entry)
{
//  AliInfo("");
  
  delete fRawEvent;
  delete fRawReader;
    
  fRawEvent = 0;
  fTree->SetBranchAddress("rawevent",&fRawEvent);
  fTree->GetEntry(entry);

  fRawReader = new AliRawReaderRoot(fRawEvent);
  
//  if ( fRawReader->GetEventHeader() )
//  {
//    fRawReader->GetEventHeader()->Print();
//  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::Notify(const char *path)
{
  AliInfo(Form("path=%s",path));
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::FinishEvent()
{
  AliInfo("");
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::Terminate()
{
  AliInfo("");
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t       
AliRawInputHandler::TerminateIO()
{
  AliInfo("");
  return kTRUE;
}
