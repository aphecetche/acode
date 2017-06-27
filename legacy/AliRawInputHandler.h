#ifndef ALIRAWINPUTHANDLER_H
#define ALIRAWINPUTHANDLER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup
/// \class AliRawInputHandler
/// \brief
/// 
/// \author Laurent Aphecetche

#ifndef ALIVEVENTHANDLER_H
#  include "AliVEventHandler.h"
#endif

class AliRawReader;
class AliRawVEvent;

class AliRawInputHandler : public AliVEventHandler
{
public:
  
  AliRawInputHandler();
  AliRawInputHandler(const char* name, const char* title);
  virtual ~AliRawInputHandler();

  virtual TTree* GetTree() const { return fTree; }
  virtual Option_t* GetDataType() const { return "RAW"; }

  virtual void         SetOutputFileName(const char* /*fname*/) {}
  virtual const char*  GetOutputFileName() { return ""; }

  virtual void         SetInputTree(TTree* tree) { fTree = tree; }

  virtual AliRawReader* GetRawReader() { return fRawReader; }
  
  // Steering 
  virtual Bool_t       Init(Option_t* opt);  
  virtual Bool_t       Init(TTree* tree, Option_t* opt);
  virtual Bool_t       BeginEvent(Long64_t entry);
  
  using AliVEventHandler::Notify;
  
  virtual Bool_t       Notify(const char *path);
  virtual Bool_t       FinishEvent();
  virtual Bool_t       Terminate();
  virtual Bool_t       TerminateIO();
  
private:
  AliRawInputHandler(const AliRawInputHandler& handler);             
  AliRawInputHandler& operator=(const AliRawInputHandler& handler);  

private:
  TTree* fTree;    //! Pointer to the tree
  AliRawReader* fRawReader; //! pointer to the rawreader
  AliRawVEvent* fRawEvent; //! pointer to rawevent
  
  ClassDef(AliRawInputHandler,1) // Implementation of AliVEventHandler for raw data
};

#endif
