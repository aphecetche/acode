#include "TGrid.h"
#include "rootFileSize.h"
#include <string>
#include <fstream>
#include "Riostream.h"
#include <map>
#include "TProof.h"
#include "TUrl.h"
#include "TFileInfo.h"
#include "TFileCollection.h"
#include "THashList.h"

const char* connect = "laphecet@nansafmaster.in2p3.fr";

int main(int argc, const char** argv)
{
  if ( argc < 3 ) 
  {
      std::cout << "usage " << argv[0] << " filename showbranches" << std::endl;
    return -1;
  }
  TString file(argv[1]);
  
  Int_t showbranches(atoi(argv[2]));
  
  if (file.Contains(".root") )
  {
    if ( file.Contains("alien://") && !gGrid ) 
    {
      TGrid::Connect("alien://");
      if (!gGrid)
      {
          std::cerr << "cannot connect to the grid" << std::endl;
        return -2;
      }
    }
    
    rootFileSize(file.Data(),(showbranches!=0));
  }
  return 0;
}
