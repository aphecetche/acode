#include "Riostream.h"
#include "TGrid.h"
#include "TFile.h"
#include "TTree.h"

void countEvents(const char* sfilelist)
{
  TGrid::Connect("alien://");
  
  ifstream filelist(sfilelist);
  char line[1024];
  Long64_t n(0);
  
  while ( filelist.getline(line,1024,'\n') )
  {
    TFile* f = TFile::Open(line);

	TTree* t = static_cast<TTree*>(f->Get("aodTree"));
	    
	n += t->GetEntries();
	
    std::cout << line << " " << t->GetEntries() << std::endl;
    
    delete f;
  }
  
  std::cout << n << " events" << std::endl;
}

int main(int argc, const char** argv)
{
  countEvents(argv[1]);
  return 0;
}
