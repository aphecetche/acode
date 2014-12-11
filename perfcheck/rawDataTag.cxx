#include "TFile.h"
#include "AliRawEventTag.h"
#include "TGrid.h"
#include "TTree.h"
#include "Riostream.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawEventHeaderVersions.h"
#include "TTimeStamp.h"
#include <map>

int main(int argc, const char** argv)
{
  if ( argc < 3 )
  {
    std::cout << "usage " << argv[0] << " raw_data_tag_file_name ref_time_stamp nseconds"<< std::endl;
    return -1;
  }
  
  TString file(argv[1]);
  
  if ( file.BeginsWith("alien://"))
  {
    TGrid::Connect("alien://");
  }
  
  TFile* f = TFile::Open(file);
  
  TString tstp(argv[2]);
  
  TTimeStamp refTimeStamp(tstp.Atoll());

  time_t timeRange(TString(argv[3]).Atoi());
  
  std::cout << refTimeStamp.AsString() << std::endl;
  
  TTree* t = static_cast<TTree*>(f->Get("T"));
  AliRawEventTag* tag(0);
  
  t->SetBranchAddress("TAG",&tag);
  
  AliRawEventHeaderV3_9::Class()->IgnoreTObjectStreamer();
	AliRawEventHeaderV3_11::Class()->IgnoreTObjectStreamer();
	AliRawEventHeaderV3_12::Class()->IgnoreTObjectStreamer();
	AliRawEventHeaderV3_13::Class()->IgnoreTObjectStreamer();
  
  Long64_t n = t->GetEntries();
  
  std::map<std::string,Long64_t> chunks;
  
  for ( Long64_t i = 0; i < n; ++i )
  {
    t->GetEntry(i);
    
    AliRawEventHeaderBase* header = tag->GetHeader();
    
    TTimeStamp ts(header->Get("Timestamp"));
    
    time_t diff = ts.GetSec() - refTimeStamp.GetSec();
    
    if ( TMath::Abs(diff) < timeRange )
    {
//      std::cout << Form("%20lld %20lld %10d %s %s",i,matchingTime,tag->GetEventNumber(),ts.AsString(),tag->GetGUID()) << std::endl;
      chunks[tag->GetGUID()]++;
    }
  }
  
  delete f;

  std::map<std::string,Long64_t>::const_iterator it;
  Long64_t nevents(0);
  
  for ( it = chunks.begin(); it != chunks.end(); ++it )
  {
    ++nevents;
    std::cout << Form("%30s %4lld",it->first.c_str(),it->second) << std::endl;
  }
  
  std::cout << Form("%20lld events in the run",n) << std::endl;
  
  std::cout << Form("Number of events in time range (%ld s): %lld",timeRange,nevents) << std::endl;

  std::cout << Form("Number of chunks having events in time range %ld",chunks.size()) << std::endl;

  return 0;
}