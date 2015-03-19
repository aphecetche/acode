void OccupancyAscii2OCDB(Int_t runNumber=67138)
{
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

  AliTestShuttle* shuttle = new AliTestShuttle(runNumber, 0, 1);
  
  const char* inputCDB = "local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB";
  
  AliTestShuttle::SetMainCDB(inputCDB);
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestReference");

  AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/OccupancyMap","alien://folder=/alice/cern.ch/user/l/laphecet/OCDB");
  
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","OCCUPANCY","MON",Form("occupancy.%d",runNumber));

  shuttle->SetInputRunType("PHYSICS");
  
  shuttle->SetDCSInput(new TMap);
  
  new AliMUONTrackerPreprocessor(shuttle);
  
  shuttle->Process();
}

