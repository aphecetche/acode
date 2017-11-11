#ifndef __COPYCDB__H__
#define __COPYCDB__H__

// Copy a OCDB objects from one CDB storage to another.

#include <vector>

void CopyCDB(Int_t runnr,
                  const char* fromURI="alien://folder=/alice/data/2017/OCDB",
                  const char* toURI="local:///alice/data/2017/OCDB/");

void CopyCDB(const char* runlist,
                  const char* fromURI="alien://folder=/alice/data/2017/OCDB",
                  const char* toURI="local:///alice/data/2017/OCDB/");

void CopyCDB(const std::vector<int>& runs,
                  const char* fromURI="alien://folder=/alice/data/2017/OCDB",
                  const char* toURI="local:///alice/data/2017/OCDB/");

// Copy and patch the MUON/Calib/Config object from one CDB storage to another.

void ChangeConfig(Int_t startRun=233721,
	              const char* fromURI="alien://folder=/alice/data/2017/OCDB",
                  const char* toURI="local:///alice/data/2017/OCDB/");

// Copy and patch the MUON/Calib/RecoParam object from one CDB storage to another.

void ChangeRecoParam(Int_t startRun=233721,
	              const char* fromURI="alien://folder=/alice/data/2017/OCDB",
                  const char* toURI="local:///alice/data/2017/OCDB/");

#endif

