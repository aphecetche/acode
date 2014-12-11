#include "dumpMC.h"

int main(int argc, char** argv)
{
  dumpMC worker(argc,argv);
  
  return worker.Exec();
}
