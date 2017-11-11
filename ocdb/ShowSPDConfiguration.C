#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliAnalysisMuMuGraphUtil.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliESDVertex.h"
#include "AliGeomManager.h"
#include "AliITSAlignMille2Module.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "Riostream.h"
#include "TBox.h"
#include "TGraph.h"
#include "TPaveText.h"
#include <algorithm>
#include <set>
#include <sstream>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TGrid.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#include <TView.h>
#include <vector>

#endif

void Draw3D(AliITSOnlineCalibrationSPDhandler *h);

struct SPDPeriod
{
  SPDPeriod() { modules.resize(240); }

  bool operator!=(const SPDPeriod& rhs)
  {
    for ( std::vector<int>::size_type i = 0; i < 240; ++i )
    {
      if ( rhs.modules[i]==0 && modules[i] != 0 ) return true;
      if ( rhs.modules[i] !=0 && modules[i] == 0 ) return true;
    }
    return false;
  }

  bool operator==(const SPDPeriod& rhs)
  {
    return !(*this != rhs);
  }
  
  int nofBadInner() const {
    int n(0);
    for ( std::vector<int>::size_type i = 0; i < 80; ++i )
    {
      if ( modules[i] ) ++n;
    }
    return n;
  };

  int nofBadOuter() const {
    int n(0);
    for ( std::vector<int>::size_type i = 80; i < 240; ++i )
    {
      if ( modules[i] ) ++n;
    }
    return n;
  };

  TH1* Create1D(const char* title) const
  {
    TH1* h = new TH1F(title,title,modules.size(),0,modules.size()-1);
    for ( std::vector<int>::size_type i = 0; i < modules.size(); ++i )
    {
      h->SetBinContent(i+1,modules[i]);
    }
    return h;
  }

  TH2F* Create2DHisto(Int_t imin, Int_t imax, const char* prefix) const;

  TH2* Create2DInner(const char* prefix) const
  {
    TString sprefix(prefix);
    sprefix += " inner ";
    return Create2DHisto(0,80,sprefix.Data());
  }

  TH2* Create2DOuter(const char* prefix) const
  {
    TString sprefix(prefix);
    sprefix += " outer ";
    return Create2DHisto(80,240,sprefix.Data());
  }

  void GetModuleIndices(Int_t module, Int_t& ix, Int_t& iy) const;

  std::vector<int> modules;
  std::vector<int> runNumbers;
};

bool operator == (const SPDPeriod& p1, const SPDPeriod& p2)
{
  return p1 == p2;
}

bool operator != (const SPDPeriod& p1, const SPDPeriod& p2)
{
  return p1 != p2;
}

bool operator < (const SPDPeriod& rhs, const SPDPeriod& other)
{
  for ( std::vector<int>::size_type i = 0; i < 240; ++i )
  {
    if ( rhs.modules[i]==0 && other.modules[i] != 0 ) return false;
    if ( rhs.modules[i] !=0 && other.modules[i] == 0 ) return false;
  }
  return true;
}

std::ostream& operator<<(std::ostream& os, const SPDPeriod& period)
{
  os << Form("%3lu runs. %2d bad inner. %2d bad outer.",period.runNumbers.size(),
             period.nofBadInner(),period.nofBadOuter());
  for ( std::vector<int>::size_type i = 0; i < period.runNumbers.size(); ++i )
  {
    if ( i % 10 == 0 ) os << std::endl << "    ";
    os << Form("%6d ",period.runNumbers[i]);
  }
  return os;
}

TCanvas* Plot(const SPDPeriod& period, const char* prefix);

//__________________________________________________________________________________________________
void PlotSPDEvolution(const std::vector<SPDPeriod>& periods)
{
  std::vector<double> vruns;
  
  for ( std::vector<SPDPeriod>::size_type i = 0; i < periods.size(); ++i )
  {
    std::cout << "Period " << i+1 << std::endl;
    
    const SPDPeriod& p = periods[i];
    
    std::cout << p << std::endl;
    
    TCanvas* c = Plot(p,Form("Period %lu / %lu - %lu runs",i+1,periods.size(),p.runNumbers.size()));
    
    c->SaveAs(Form("period-%lu.pdf",i+1));
    
    for ( std::vector<int>::size_type j = 0; j < p.runNumbers.size(); ++j )
    {
      vruns.push_back(p.runNumbers[j]);
    }
  }
  
  std::sort(vruns.begin(),vruns.end());
  
  std::cout << "---- Total "  << vruns.size() << " runs grouped into " << periods.size() << " periods" << std::endl;
  
}

//______________________________________________________________________________
void ReadIntegers(const char* filename, std::vector<int>& integers, Bool_t resetVector)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  /// copied from AliAnalysisTriggerScalers::ReadIntegers (it's a bad thing to copy&paste,
  /// but I did not want to introduce more dependencies to this macro)
  
  if ( gSystem->AccessPathName(gSystem->ExpandPathName(filename))==kTRUE )
  {
    return;
  }
  std::ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
  std::set<int> runset;
  
  if (!resetVector)
  {
    for ( std::vector<int>::size_type j = 0; j < integers.size(); ++ j )
    {
      runset.insert(integers[j]);
    }
  }
  
  char line[10000];
  
  in.getline(line,10000,'\n');
  
  TString sline(line);
  
  if (sline.Contains(","))
  {
    TObjArray* a = sline.Tokenize(",");
    TIter next(a);
    TObjString* s;
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }
    delete a;
  }
  else
  {
    runset.insert(sline.Atoi());
    
    while ( in >> i )
    {
      runset.insert(i);
    }
  }
  
  integers.clear();
  
  for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
  {
    integers.push_back((*it));
  }
  
  std::sort(integers.begin(),integers.end());
}

//______________________________________________________________________________
TH2F* SPDPeriod::Create2DHisto(Int_t imin, Int_t imax, const char* prefix) const
{
  std::vector<Double_t> x;
  std::vector<Double_t> y;

  TGeoHMatrix matrix;
  
  Double_t localBottomLeft[3];
  Double_t localTopRight[3];
  
  if ( imin == 0 )
  {
    localBottomLeft[0] = -0.6375;
    localBottomLeft[1] = 0.0;
    localBottomLeft[2] = -3.48;
    localTopRight[0] = 0.6375;
    localTopRight[1] = 0.0;
    localTopRight[2] = 3.48;
  }
  else
  {
    localBottomLeft[0] = -0.6375;
    localBottomLeft[1] = 0.0;
    localBottomLeft[2] = 3.48;
    localTopRight[0] = 0.6375;
    localTopRight[1] = 0.0;
    localTopRight[2] = -3.48;
  }
  
  Double_t globalBottomLeft[3];
  Double_t globalTopRight[3];
  
  int nx=4;
  
  // first get x-axis (representing module's z) binning
  
  for (Int_t i=imin+nx-1; i>=imin; --i)
  {
    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(i);
    AliITSAlignMille2Module::SensVolMatrix(vid,&matrix);
    
    memset(globalBottomLeft,0,sizeof(Double_t)*3);
    
    matrix.LocalToMaster(localBottomLeft,globalBottomLeft);
    matrix.LocalToMaster(localTopRight,globalTopRight);

//    std::cout << Form("module %i xl = %e xr = %e",i,globalBottomLeft[2],globalTopRight[2]) << std::endl;
    

    x.push_back(globalBottomLeft[2]);

  }

  // must get the upper edge for the last bin

  memset(globalTopRight,0,sizeof(Double_t)*3);
  matrix.LocalToMaster(localTopRight,globalTopRight);
  x.push_back(globalTopRight[2]);

  // now get y-axis (representing module's phi) binning
  
  Double_t phiDown(0);
  
  for (Int_t i=imin; i<imax; i+=nx)
  {
    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(i);
    AliITSAlignMille2Module::SensVolMatrix(vid,&matrix);
    
    memset(globalBottomLeft,0,sizeof(Double_t)*3);
    
    matrix.LocalToMaster(localBottomLeft,globalBottomLeft);
    phiDown = TMath::ATan2(globalBottomLeft[1],globalBottomLeft[0]);
    if ( phiDown < 0 )
    {
      phiDown += 2*TMath::Pi();
    }

    y.push_back(phiDown);
    
  }

  memset(globalTopRight,0,sizeof(Double_t)*3);
  matrix.LocalToMaster(localTopRight,globalTopRight);
  Double_t phiUp = TMath::ATan2(globalTopRight[1],globalTopRight[0]);
  if ( phiUp < 0 )
  {
    phiUp += 2*TMath::Pi();
  }
  if ( phiUp < phiDown )
  {
    phiUp  = 2*TMath::Pi(); // dirty shortcut to get the phimax fit in the histogram...
  }
  
  y.push_back(phiUp);

  TString sname(prefix);
  
  sname += " 2D map of inactive modules ";
  TH2F* h = new TH2F(sname.Data(),sname.Data(),x.size()-1,&x[0],y.size()-1,&y[0]);
  
  h->GetXaxis()->SetTickLength(0.0);
  h->GetYaxis()->SetTickLength(0.0);
  
  int ix,iy;

  for (Int_t module=imin; module<imax; module++)
  {
    if ( modules[module] )
    {
      GetModuleIndices(module,ix,iy);
      h->SetBinContent(ix,iy,h->GetBinContent(ix,iy)+1.0);
    }
  }

  return h;
}

//______________________________________________________________________________
void ShowSPDConfiguration(Int_t runNb=0, const char *ocdblocation="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2012/OCDB", bool grid=kFALSE,bool threed=kFALSE){


  gStyle->SetOptStat(0);

  if(grid){
    TGrid::Connect("alien://");
    if(!gGrid){
      printf("no grid connection is available, exiting.\n");
      return;
    }
  }
  AliCDBManager::Instance();
  AliCDBManager::Instance()->SetDefaultStorage(ocdblocation);
  AliCDBManager::Instance()->SetRun(runNb);
  AliGeomManager::LoadGeometry();

  AliITSOnlineCalibrationSPDhandler *h = new AliITSOnlineCalibrationSPDhandler();
  h->ReadDeadFromDB(runNb,ocdblocation);

  if(threed) {
    Draw3D(h);  
    return;
  }

  Double_t nact[2]={0.,0.};

  TCanvas *c = new TCanvas("c","Active Modules ",500,700);
  c->Divide(1,2);
  c->cd(1);
  TH2D *hPhiZInner = new TH2D("hPhiZInner","Inner layer Active Modules ",200,-20,20,3,0,2*TMath::Pi());
  hPhiZInner->SetXTitle("Z (cm)");
  hPhiZInner->SetYTitle("#varphi (rad)");
  hPhiZInner->Draw();

  for(Int_t i=0; i<80; ++i){
    if(h->GetNrBad(i) <1) {
      TGeoHMatrix matrix;
      int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(i);
      AliITSAlignMille2Module::SensVolMatrix(vid,&matrix);
      Double_t local0[3],local1[3],local2[3],local3[3]; // local position of the four angles
      local0[0]=-0.6375; local0[1]=0; local0[2]= 3.48; 
      local1[0]=-0.6375; local1[1]=0; local1[2]= -3.48; 
      local2[0]=0.6375; local2[1]=0; local2[2]= -3.48; 
      local2[0]=0.6375; local3[1]=0; local3[2]= 3.48; 
      Double_t global0[3],global1[3],global2[3],global3[3];
      matrix.LocalToMaster(local0,global0);
      matrix.LocalToMaster(local1,global1);
      matrix.LocalToMaster(local2,global2);
      matrix.LocalToMaster(local3,global3);
      Double_t phiUp = atan2(global0[1],global0[0]);
      if(phiUp<0) phiUp+=2*TMath::Pi();
      Double_t phiDown = atan2(global2[1],global2[0]);
      if(phiDown<0) phiDown+=2*TMath::Pi();
      int color=kBlue;
      if ( i == 30 ) color = kRed;
      TLine *lhor1 = new TLine(global0[2],phiDown,global1[2],phiDown); lhor1->Draw("same");
      lhor1->SetLineColor(color);
      lhor1->SetLineWidth(3);
      TLine *lver1 = new TLine(global1[2],phiDown,global2[2],phiUp); lver1->Draw("same");
      lver1->SetLineColor(color);
      lver1->SetLineWidth(3); 
      TLine *lhor2 = new TLine(global2[2],phiUp,global3[2],phiUp); lhor2->Draw("same");
      lhor2->SetLineColor(color);
      lhor2->SetLineWidth(3); 
      TLine *lver2 = new TLine(global3[2],phiUp,global0[2],phiDown); lver2->Draw("same");
      lver2->SetLineColor(color);
      lver2->SetLineWidth(3);
      nact[0]++;
    }
    else
    {
      std::cout << "Inner module # " << i << " is missing" << std::endl;
    }
  }
  c->cd(2);

  TH2D *hPhiZOuter = new TH2D("hPhiZOuter","Outer layer Active Modules ",200,-20,20,3,0,2*TMath::Pi());
  hPhiZOuter->SetXTitle("Z (cm)");
  hPhiZOuter->SetYTitle("#varphi (rad)");
  hPhiZOuter->Draw();
  for(Int_t i=80; i<240; i++){
    TGeoHMatrix matrix;
    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(i);
    AliITSAlignMille2Module::SensVolMatrix(vid,&matrix);
    Double_t local[4][3]; // local position of the four angles
    local[0][0]=0.6375; local[0][1]=0; local[0][2]= -3.48;
    local[1][0]=0.6375; local[1][1]=0; local[1][2]= 3.48;
    local[2][0]=-0.6375; local[2][1]=0; local[2][2]= 3.48;
    local[2][0]=-0.6375; local[3][1]=0; local[3][2]= -3.48;
    Double_t global[4][3];
    for(Int_t j=0; j<4; j++){
      matrix.LocalToMaster(local[j],global[j]);
    }
    Double_t phiUp = atan2(global[0][1],global[0][0]);
    if(phiUp<0) phiUp+=2*TMath::Pi();
    Double_t phiDown = atan2(global[2][1],global[2][0]);
    if(phiDown<0) phiDown+=2*TMath::Pi();
    if(i>235) if(phiUp < 0.1) phiUp = TMath::Pi()*2;
    //  printf("module %i  -   phiDown %f   phiUp %f \n",i,phiDown,phiUp); 
    if((h->GetNrBad(i))<1) {
      TLine *lhor1 = new TLine(global[0][2],phiUp,global[1][2],phiUp); lhor1->Draw("same");
      lhor1->SetLineColor(kBlue);
      lhor1->SetLineWidth(2); 
      TLine *lver1 = new TLine(global[1][2],phiUp,global[2][2],phiDown); lver1->Draw("same");
      lver1->SetLineColor(kBlue);
      lver1->SetLineWidth(2); 
      TLine *lhor2 = new TLine(global[2][2],phiDown,global[3][2],phiDown); lhor2->Draw("same");
      lhor2->SetLineColor(kBlue);
      lhor2->SetLineWidth(2); 
      TLine *lver2 = new TLine(global[3][2],phiDown,global[0][2],phiUp); lver2->Draw("same");
      lver2->SetLineColor(kBlue);
      lver2->SetLineWidth(2); 
      nact[1]++; 
    } 
  }
  printf("  \n   Number of Active SPD modules (->Total)  : inner %3.0f (80) %f   -  outer %3.0f (160)  %f \n",nact[0],nact[0]/80.,nact[1],nact[1]/160.);
  c->SaveAs(Form("active%i.png",runNb));
}


//______________________________________________________________________________
void Draw3D(AliITSOnlineCalibrationSPDhandler *h){

  TGeoHMatrix m2t[240];
  for(Int_t imod=0; imod<240; imod++){
    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(imod);
    AliITSAlignMille2Module::SensVolMatrix(vid,&m2t[imod]);
  }

  delete gGeoManager;

  new TGeoManager("SPD","active");

  TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
  TGeoMedium *none = new TGeoMedium("Vacuum",0,vacuum);
  TGeoVolume *top = gGeoManager->MakeBox("TOP",none,500,500,500);
  gGeoManager->SetTopVolume(top);

  TGeoVolume *ladder = gGeoManager->MakeBox("ladder",none,0.6375,0.001/2,3.48);

  Int_t nActive[2]={0,0};
  for(Int_t imod=0; imod<240; imod++){
    TGeoRotation *rot  = new TGeoRotation();
    rot->SetMatrix(m2t[imod].GetRotationMatrix());
    TGeoCombiTrans *matrix = new TGeoCombiTrans(m2t[imod].GetTranslation()[0],m2t[imod].GetTranslation()[1],m2t[imod].GetTranslation()[2],rot);
    if((40960-h->GetNrBad(imod))>0) {
      top->AddNode(ladder,imod,matrix);
      if(imod<80) nActive[0]++;
      else nActive[1]++;
    }
  }

  printf("  \n\n   Number of Active SPD modules (->Total)  : inner %i (80) outer %i (160) \n\n\n",nActive[0],nActive[1]);
  gGeoManager->CloseGeometry();
  top->Draw("ogl");
  gPad->GetView()->ShowAxis();

}

//______________________________________________________________________________
void SPDPeriod::GetModuleIndices(Int_t module, Int_t& ix, Int_t& iy) const
{
  Int_t firstModule = 0;
  
  if (module>=80) firstModule=80;
  
  iy = (module-firstModule)/4;
  ix = (module-firstModule) - iy*4;
  
  ++iy;
  ++ix;
}

//______________________________________________________________________________
TCanvas* Plot(const SPDPeriod& period, const char* prefix)
{
  TCanvas* c = new TCanvas(prefix,prefix);
  c->Divide(1,2);
  c->cd(1);
  
  TVirtualPad* pad = gPad;
  gPad->Divide(2,1);
  
  pad->cd(1);
  
  TH2* hinner = period.Create2DInner(prefix);
  hinner->Draw("colz");
  
  pad->cd(2);

  TH2* houter = period.Create2DOuter(prefix);
  houter->Draw("colz");

  c->cd(2);
  pad = gPad;
  gPad->Divide(2,1);
  pad->cd(1);
  
  TH1* h = period.Create1D(prefix);
  
  h->Draw("hist");
  
  pad->cd(2);
  
  gPad->Draw();
  
  TPaveText* text = new TPaveText(0,0,1,1,"NDC");
  text->SetFillStyle(0);
  text->Draw();
  text->SetTextSize(0.03);
  std::ostringstream* runs = new std::ostringstream;
  int n(0);
  
  for ( std::vector<int>::size_type i = 0; i < period.runNumbers.size(); ++i )
  {
    (*runs) << period.runNumbers[i] << ",";
    ++n;
    if ( n == 10 )
    {
      text->AddText(runs->str().c_str());
      runs = new std::ostringstream;
      n = 0;
    }
  }
  
  if ( runs->str().length() > 0 )
  {
    text->AddText(runs->str().c_str());
  }
  
  c->Draw();
  
  return c;
}

//______________________________________________________________________________
Bool_t SPDConfigEvolution(const char* runlist, const char* ocdbPath)
{
  gStyle->SetOptStat(0);

  std::vector<int> runNumbers;
  
  ReadIntegers(runlist,runNumbers,kTRUE);
  
  if (runNumbers.empty())
  {
    std::cerr << "Cannot work with an empty run list !" << std::endl;
    return kFALSE;
  }
  
  std::vector<SPDPeriod> periods;
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath);
  man->SetRun(runNumbers[0]); // here we assume the geometry is not changing
  // within the period described by the runlist...
  AliGeomManager::LoadGeometry();
  
  std::vector<double> vruns;
  std::vector<double> vbadsInner;
  std::vector<double> vbadsOuter;

  for ( std::vector<int>::size_type i = 0; i < runNumbers.size(); ++i )
  {
    Int_t runNumber = runNumbers[i];
    
    vruns.push_back(runNumber);

    SPDPeriod thisPeriod;
    
    thisPeriod.runNumbers.push_back(runNumber);
    
    AliCDBEntry* cdbEntry = man->Get("ITS/Calib/SPDDead", runNumber);
    
    if (!cdbEntry)
    {
      Warning("AliITSOnlineCalibrationSPDhandler::ReadDeadFromDB","Calibration for run %d not found in database.",runNumber);
      continue;
    }

    TObjArray* spdEntry = static_cast<TObjArray*>(cdbEntry->GetObject());
    if (spdEntry)
    {
      for (Int_t module=0; module<240; ++module)
      {
        AliITSCalibrationSPD* calibSPD = static_cast<AliITSCalibrationSPD*>(spdEntry->At(module));
        
        if (!calibSPD)
        {
          std::cerr << "got a null calibSPD for module " << module << " ! " << std::endl;
          continue;
        }
        
        if ( calibSPD->GetNrBad() )
        {
          thisPeriod.modules[module]++;
        }
      }
      
      std::vector<SPDPeriod>::iterator it = std::find(periods.begin(),periods.end(),thisPeriod);//periods.find(thisPeriod);
      
      if ( it != periods.end() )
      {
        it->runNumbers.push_back(runNumber);
        for (Int_t module=0; module<240; ++module)
        {
          it->modules[module] += thisPeriod.modules[module];
        }
      }
      else
      {
        periods.push_back(thisPeriod);
      }

      vbadsInner.push_back(thisPeriod.nofBadInner());
      vbadsOuter.push_back(thisPeriod.nofBadOuter());
    }
  }
  
  PlotSPDEvolution(periods);
  
  TCanvas* canvas = new TCanvas("Number of inactive SPD modules","Number of inactive SPD modules");

  TGraph* ginner = new TGraph(vruns.size(),&vruns[0],&vbadsInner[0]);
  TGraph* gouter = new TGraph(vruns.size(),&vruns[0],&vbadsOuter[0]);
 
   if (!ginner || !gouter) return kFALSE;
 
   ginner->SetName("Inner Inactive Modules");
  gouter->SetName("Outer Inactive Modules");
 
 
  AliAnalysisMuMuGraphUtil gu;
  gu.ShouldDrawPeriods(kTRUE);
  TObjArray graphs;
  graphs.Add(ginner);
  graphs.Add(gouter);
 
  Double_t ymin,ymaxin,ymaxout;
 
  gu.GetYMinAndMax(*ginner, ymin, ymaxin);
  gu.GetYMinAndMax(*gouter, ymin, ymaxout);
 
  gu.PlotSameWithLegend(graphs,0,TMath::Max(ymaxin,ymaxout)+10);

  canvas->Print("inactive.spd.pdf");
  
  return kTRUE;
}

//______________________________________________________________________________
void MeanVertexEvolution(const char* runlist, const char* ocdbPath)
{
  gStyle->SetOptStat(0);
  
  std::vector<int> runNumbers;
  
  ReadIntegers(runlist,runNumbers,kTRUE);
  
  if (runNumbers.empty())
  {
    std::cerr << "Cannot work with an empty run list !" << std::endl;
    return;
  }
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath);
  man->SetRun(runNumbers[0]);

  
  std::vector<double> vruns;
  std::vector<double> vmean;
  std::vector<double> vsigma;
  
  for (std::vector<int>::size_type i = 0; i < runNumbers.size(); ++i )
  {
    int run = runNumbers[i];
    AliCDBEntry* entry = man->Get("GRP/Calib/MeanVertexSPD",run);
    if (!entry) continue;
    AliESDVertex* v = static_cast<AliESDVertex*>(entry->GetObject());
    vruns.push_back(run);
    vmean.push_back(v->GetZ()*10);
    vsigma.push_back(v->GetZRes());
  }
  
  TGraph* gmean = new TGraph(vruns.size(),&vruns[0],&vmean[0]);
  gmean->SetName("Mean vertex (mm)");
  TGraph* gsigma = new TGraph(vruns.size(),&vruns[0],&vsigma[0]);
  gsigma->SetName("Sigma vertex (cm)");

  AliAnalysisMuMuGraphUtil gu;
  gu.ShouldDrawPeriods(kTRUE);
  TObjArray graphs;
  graphs.Add(gmean);
  graphs.Add(gsigma);
  
  gmean->Print();
  
//  gmean->GetXaxis()->SetNoExponent();
//  gsigma->GetXaxis()->SetNoExponent();
  
//  TCanvas* c = gu.DrawWith2Scales(*gmean, *gsigma, "mean vertex evolution");
//  c->Draw();

  gu.PlotSameWithLegend(graphs,-10,40);
}

