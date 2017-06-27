#ifndef AliAnalysisTaskRaw_h
#define AliAnalysisTaskRaw_h

class AliMUONTrackerDataMaker;

#include "AliAnalysisTask.h"
#include "TString.h"

class AliAnalysisTaskRaw : public AliAnalysisTask {
 public:
  AliAnalysisTaskRaw(const char *name = "AliAnalysisTaskRaw", 
                     const char* ocdbpath="alien://folder=/alice/data/2009/OCDB");
  virtual ~AliAnalysisTaskRaw() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  virtual void LocalInit();
  
private:
  AliMUONTrackerDataMaker* fDataMaker;
  TString fOCDBPath;
  
  AliAnalysisTaskRaw(const AliAnalysisTaskRaw&); // not implemented
  AliAnalysisTaskRaw& operator=(const AliAnalysisTaskRaw&); // not implemented
  
  ClassDef(AliAnalysisTaskRaw, 1); // example of analysis
};

#endif
