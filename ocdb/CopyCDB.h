#ifndef __COPYCDB__H__
#define __COPYCDB__H__

// Copy a OCDB objects from one CDB storage to another.

#include <vector>

void CopyCDB(Int_t runnr,
                  const char* fromURI="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB",
                  const char* toURI="local:///alice/cern.ch/user/l/laphecet/OCDB/");

void CopyCDB(const char* runlist,
                  const char* fromURI="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB",
                  const char* toURI="local:///alice/cern.ch/user/l/laphecet/OCDB/");

void CopyCDB(const std::vector<int>& runs,
                  const char* fromURI="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB",
                  const char* toURI="local:///alice/cern.ch/user/l/laphecet/OCDB/");

// Copy and patch the MUON/Calib/Config object from one CDB storage to another.

void ChangeConfig(Int_t startRun=233721,
	              const char* fromURI="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB",
                  const char* toURI="local:///alice/cern.ch/user/l/laphecet/OCDB/");

// Copy and patch the MUON/Calib/RecoParam object from one CDB storage to another.

void ChangeRecoParam(Int_t startRun=233721,
	              const char* fromURI="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB",
                  const char* toURI="local:///alice/cern.ch/user/l/laphecet/OCDB/");

#endif

