///
/// Benchmark for Root
///
/// adapted from root v5-34-30 $ROOTSYS/src/test/stress*.cxx 
/// in order to be able to more easily run them in //
///

#include "Event.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "Riostream.h"
#include "rootstress.h"
#include "TApplication.h"
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TDecompBK.h"
#include "TDecompChol.h"
#include "TDecompLU.h"
#include "TDecompQRH.h"
#include "TDecompSVD.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TH2.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDLazy.h"
#include "TMatrixDSparse.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDUtils.h"
#include "TMatrixF.h"
#include "TMatrixFLazy.h"
#include "TMatrixFSparse.h"
#include "TMatrixFSym.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include <Compression.h>
#include <Riostream.h>
#include <stdlib.h>
#include <TApplication.h>
#include <TArrayD.h>
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClassTable.h>
#include <TCut.h>
#include <TCutG.h>
#include <TEventList.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TPostScript.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTreeCache.h>
#include <map>

Double_t gToleranceMult = 1.e-3;

#define EPSILON 1.0e-14
#define VERBOSE 0

std::string gInputFiles;

Int_t gVerbose      = 0;
Int_t gNrLoop;
const Int_t nrSize  = 20;
const Int_t gSizeA[] = {5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,300,500,700,1000};

//------------------------------------------------------------------------
void StatusPrint(Int_t id,const TString &title, Int_t nsuccess, Int_t nattempts)
{
#if VERBOSE
  // Print test program number and its title
  const Int_t kMAX = 65;
  Char_t number[4];
  snprintf(number,4,"%2d",id);
  TString header = TString("Test ")+number+" : "+title;
  const Int_t nch = header.Length();
  for (Int_t i = nch; i < kMAX; i++) header += '.';
  std::cout << header << " " << nsuccess << " out of " << nattempts << std::endl;
#endif
}

//------------------------------------------------------------------------
void StatusPrint(Int_t id,const TString &title,Bool_t status)
{
#if VERBOSE
  // Print test program number and its title
  const Int_t kMAX = 65;
  Char_t number[4];
  snprintf(number,4,"%2d",id);
  TString header = TString("Test ")+number+" : "+title;
  const Int_t nch = header.Length();
  for (Int_t i = nch; i < kMAX; i++) header += '.';
  cout << header << (status ? "OK" : "FAILED") << endl;
#endif
}

//______________________________________________________________________________
Double_t RosenBrock(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

//______________________________________________________________________________
Bool_t RunRosenBrock()
{
//
// F(x,y) = 100 (y-x^2)^2 + (1-x)^2
//
//   start point: F(-1.2,1.0) = 24.20
//   minimum    : F(1.0,1.0)  = 0.
//
// This narrow, parabolic valley is probably the best known of all test cases. The floor
// of the valley follows approximately the parabola y = x^2+1/200 .
// There is a region where the covariance matrix is not positive-definite and even a path
// where it is singular . Stepping methods tend to perform at least as well as gradient
//  method for this function .
// [Reference: Comput. J. 3,175 (1960).]

  Bool_t ok = kTRUE;
  const int nvars = 2;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&RosenBrock, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "x", -1.2, step);
  min->SetVariable(1, "y",  1.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

//______________________________________________________________________________
Double_t Wood4(const Double_t *par)
{
  const Double_t w = par[0]; 
  const Double_t x = par[1]; 
  const Double_t y = par[2]; 
  const Double_t z = par[3]; 

  const Double_t w1 = w-1;
  const Double_t x1 = x-1;
  const Double_t y1 = y-1;
  const Double_t z1 = z-1;
  const Double_t tmp1 = x-w*w;
  const Double_t tmp2 = z-y*y;

  return 100*tmp1*tmp1+w1*w1+90*tmp2*tmp2+y1*y1+10.1*(x1*x1+z1*z1)+19.8*x1*z1;
}

//______________________________________________________________________________
Bool_t RunWood4()
{
//
// F(w,x,y,z) = 100 (y-w^2)^2 + (w-1)^2 + 90 (z-y^2)^2
//              + (1-y)^2 + 10.1 [(x-1)^2 + (z-1)^2]
//              + 19.8 (x-1)(z-1)
//
//   start point: F(-3,-1,-3,-1) = 19192
//   minimum    : F(1,1,1,1)  =   0.
//
// This is a fourth-degree polynomial which is reasonably well-behaved near the minimum,
// but in order to get there one must cross a rather flat, four-dimensional "plateau"
// which often causes minimization algorithm to get "stuck" far from the minimum. As
// such it is a particularly good test of convergence criteria and simulates quite well a
// feature of many physical problems in many variables where no good starting
// approximation is known .
// [Reference: Unpublished. See IBM Technical Report No. 320-2949.]


  Bool_t ok = kTRUE;
  const int nvars = 4;
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  if (!min)
    {
      std::cerr << "RunWood4(): failed to create ROOT::Math::Minimizer" << std::endl;
      return 0;
    }
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&Wood4, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "w", -3.0, step);
  min->SetVariable(1, "x", -1.0, step);
  min->SetVariable(2, "y", -3.0, step);
  min->SetVariable(3, "z", -1.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;
  delete min;

  return ok;
}

//______________________________________________________________________________
Double_t Powell(const Double_t *par)
{
  const Double_t w = par[0]; 
  const Double_t x = par[1]; 
  const Double_t y = par[2]; 
  const Double_t z = par[3]; 

  const Double_t tmp1 = w+10*x;
  const Double_t tmp2 = y-z;
  const Double_t tmp3 = x-2*y;
  const Double_t tmp4 = w-z;

  return tmp1*tmp1+5*tmp2*tmp2+tmp3*tmp3*tmp3*tmp3+10*tmp4*tmp4*tmp4*tmp4;
}

//______________________________________________________________________________
Bool_t RunPowell()
{
//
// F(w,x,y,z) = (w+10x)^2 + 5(y-z)^2 + (x-2y)^4+ 10 (w-z)^4
//
//   start point: F(-3,-1,0,1) = 215
//   minimum    : F(0,0,0,0)  =   0.
//
// This function is difficult because its matrix of second derivatives becomes singular
//  at the minimum. Near the minimum the function is given by (w + 10x)^2 + 5 (y-5)^2
// which does not determine the minimum uniquely.
// [Reference: Comput. J. 5, 147 (1962).]

  Bool_t ok = kTRUE;
  const int nvars = 4;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&Powell, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "w", +3.0, step);
  min->SetVariable(1, "x", -1.0, step);
  min->SetVariable(2, "y",  0.0, step);
  min->SetVariable(3, "z", +1.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

//______________________________________________________________________________
Double_t Fletcher(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];
  const Double_t z = par[2];

  Double_t psi;
  if (x > 0)
    psi = TMath::ATan(y/x)/2/TMath::Pi();
  else if (x < 0)
    psi = 0.5+TMath::ATan(y/x)/2/TMath::Pi();
  else
    psi = 0.0;

  const Double_t tmp1 = z-10*psi;
  const Double_t tmp2 = TMath::Sqrt(x*x+y*y)-1;

  return 100*(tmp1*tmp1+tmp2*tmp2)+z*z;
}

//______________________________________________________________________________
Bool_t RunFletcher()
{
//
// F(x,y,z) = 100 {[z - 10 G(x,y)]^2 + ( (x^2+y^2)^1/2 - 1 )^2} + z^2
//
//                     | arctan(y/x)        for x > 0
// where 2 pi G(x,y) = |
//                     | pi + arctan(y/x)   for x < 0
//
//   start point: F(-1,0,0) = 2500
//   minimum    : F(1,0,0)  =   0.
//
// F is defined only for -0.25 < G(x,y) < 0.75
//
// This is a curved valley problem, similar to Rosenbrock's, but in three dimensions .
// [Reference: Comput. J. 6, 163 (1963).]

  Bool_t ok = kTRUE;
  const int nvars = 3;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&Fletcher, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "x", -1.0, step);
  min->SetVariable(1, "y",  0.0, step);
  min->SetVariable(2, "z",  0.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

//______________________________________________________________________________
Double_t GoldStein1(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];

  const Double_t tmp1 = x+y+1;
  const Double_t tmp2 = 19-14*x+3*x*x-14*y+6*x*y+3*y*y;
  const Double_t tmp3 = 2*x-3*y;
  const Double_t tmp4 = 18-32*x+12*x*x+48*y-36*x*y+27*y*y;

  return (1+tmp1*tmp1*tmp2)*(30+tmp3*tmp3*tmp4);
}

//______________________________________________________________________________
Bool_t RunGoldStein1()
{
//
// F(x,y) = (1 + (x+y+1)^2 * (19-14x+3x^2-14y+6xy+3y^2))
//           * (30 + (2x-3y)^2 * (18-32x+12x^2+48y-36xy+27y^2))
// 
//   start point     : F(-0.4,-0,6) = 35
//   local  minima   : F(1.2,0.8)   = 840
//                     F(1.8,0.2)   = 84
//                     F(-0.6,-0.4) = 30
//   global minimum  : F(0.0,-1.0)  = 3
//   
// This is an eighth-order polynomial in two variables which is well behaved near each
// minimum, but has four local minima and is of course non-positive-definite in many
// regions. The saddle point between the two lowest minima occurs at F(-0.4,-0.6)=35
// making this an interesting start point .
// [Reference: Math. Comp. 25, 571 (1971).]

  Bool_t ok = kTRUE;
  const int nvars = 2;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&GoldStein1, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -2.0;
  double xmax =  2.0;
  double ymin = 3.0;
  double tolerance = 1.e-3;

  min->SetLimitedVariable(0, "x", -0.3999, step, xmin, xmax);
  min->SetLimitedVariable(1, "y", -0.6,    step, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > ymin + tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

//______________________________________________________________________________
Double_t GoldStein2(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];

  const Double_t tmp1 = x*x+y*y-25;
  const Double_t tmp2 = TMath::Sin(4*x-3*y);
  const Double_t tmp3 = 2*x+y-10;

  return TMath::Exp(0.5*tmp1*tmp1)+tmp2*tmp2*tmp2*tmp2+0.5*tmp3*tmp3;
}

//______________________________________________________________________________
Bool_t RunGoldStein2()
{
//
// F(x,y) = (1 + (x+y+1)^2 * (19-14x+3x^2-14y+6xy+3y^2))
//           * (30 + (2x-3y)^2 * (18-32x+12x^2+48y-36xy+27y^2))
// 
//   start point     : F(1.6,3.4) =
//   global minimum  : F(3,4)     = 1
//   
// This function has many local minima .
// [Reference: Math. Comp. 25, 571 (1971).]

  Bool_t ok = kTRUE;
  const int nvars = 2;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&GoldStein2, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double ymin = 1.0;
  double tolerance = 1.e-2;

  min->SetLimitedVariable(0, "x", +1.0, step, xmin, xmax);
  min->SetLimitedVariable(1, "y", +3.2, step, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > ymin + tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

Double_t seed = 3;
Int_t  nf;
TMatrixD A;
TMatrixD B;
TVectorD x0;
TVectorD sx0;
TVectorD cx0;
TVectorD sx;
TVectorD cx;
TVectorD v0;
TVectorD v;
TVectorD r;

//______________________________________________________________________________
Double_t TrigoFletcher(const Double_t *par)
{
  Int_t i;
  for (i = 0; i < nf ; i++) {
    cx0[i] = TMath::Cos(x0[i]);
    sx0[i] = TMath::Sin(x0[i]);
    cx [i] = TMath::Cos(par[i]);
    sx [i] = TMath::Sin(par[i]);
  }

  v0 = A*sx0+B*cx0;
  v  = A*sx +B*cx;
  r  = v0-v;
 
  return r * r;
}

//______________________________________________________________________________
Bool_t RunTrigoFletcher()
{
//
// F(\vec{x}) = \sum_{i=1}^n ( E_i - \sum_{j=1}^n (A_{ij} \sin x_j + B_{ij} \cos x_j) )^2
// 
//   where E_i = \sum_{j=1}^n ( A_{ij} \sin x_{0j} + B_{ij} \cos x_{0j} )
//
//   B_{ij} and A_{ij} are random matrices composed of integers between -100 and 100;
//   for j = 1,...,n: x_{0j} are any random numbers, -\pi < x_{0j} < \pi;
//
//   start point : x_j = x_{0j} + 0.1 \delta_j,  -\pi < \delta_j < \pi
//   minimum     : F(\vec{x} = \vec{x}_0) = 0
//   
// This is a set of functions of any number of variables n, where the minimum is always
// known in advance, but where the problem can be changed by choosing different
// (random) values of the constants A_{ij}, B_{ij}, and x_{0j} . The difficulty can be
// varied by choosing larger starting deviations \delta_j . In practice, most methods
// find the "right" minimum, corresponding to \vec{x} = \vec{x}_0, but there are usually
// many subsidiary minima.
// [Reference: Comput. J. 6 163 (1963).]

 
  const Double_t pi = TMath::Pi();
  Bool_t ok = kTRUE;
  Double_t delta = 0.1;
  Double_t tolerance = 1.e-2;

  for (nf = 5; nf<32;nf +=5) {
     ROOT::Math::Minimizer* min =
       ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
     min->SetPrintLevel(gVerbose);

     ROOT::Math::Functor f(&TrigoFletcher, nf);
     min->SetFunction(f);

     A.ResizeTo(nf,nf);
     B.ResizeTo(nf,nf);
     x0.ResizeTo(nf);
     sx0.ResizeTo(nf);
     cx0.ResizeTo(nf);
     sx.ResizeTo(nf);
     cx.ResizeTo(nf);
     v0.ResizeTo(nf);
     v.ResizeTo(nf);
     r.ResizeTo(nf);
     A.Randomize(-100.,100,seed);
     B.Randomize(-100.,100,seed);
     for (Int_t i = 0; i < nf; i++) {
       for (Int_t j = 0; j < nf; j++) {
         A(i,j) = Int_t(A(i,j));
         B(i,j) = Int_t(B(i,j));
       }
     }

     x0.Randomize(-pi,pi,seed);
     TVectorD x1(nf); x1.Randomize(-delta*pi,delta*pi,seed);
     x1+= x0;

     for (Int_t i = 0; i < nf; i++)
       min->SetLimitedVariable(i, Form("x_%d",i), x1[i], 0.01, -pi*(1+delta), +pi*(1+delta));

     min->SetMaxFunctionCalls(100000);
     min->SetTolerance(tolerance * gToleranceMult);

     min->Minimize();

     if (min->MinValue() > tolerance)
       ok = kFALSE;

     delete min;
  }  

  return ok;
}

//______________________________________________________________________________
Double_t stressFit(const char *type, const char *algo, Int_t N)
{
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(type, algo);

#if VERBOSE
  std::cout << "******************************************************************" <<std::endl;
  std::cout << "*  Minimization - S T R E S S suite                              *" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;
#endif

   TStopwatch timer;
   timer.Start();

#if VERBOSE
  std::cout << "*  Starting  S T R E S S  with fitter : "
            << ROOT::Math::MinimizerOptions::DefaultMinimizerType() << " / "
            << ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo() << std::endl;
  std::cout << "******************************************************************" << std::endl;
#endif

  gBenchmark->Start("stressFit");

  int okRosenBrock    = 0;
  int okWood          = 0;
  int okPowell        = 0;
  int okFletcher      = 0;
  int okGoldStein1    = 0;
  int okGoldStein2    = 0;
  int okTrigoFletcher = 0;
  Int_t i;
  for (i = 0; i < N; i++) if (RunWood4()) okWood++;
  StatusPrint(1, "Wood", okWood, N);
  for (i = 0; i < N; i++) if (RunRosenBrock()) okRosenBrock++;
  StatusPrint(2, "RosenBrock", okRosenBrock, N);
  for (i = 0; i < N; i++) if (RunPowell()) okPowell++;
  StatusPrint(3, "Powell", okPowell, N);
  for (i = 0; i < N; i++) if (RunFletcher()) okFletcher++;
  StatusPrint(4, "Fletcher", okFletcher, N);
  for (i = 0; i < N; i++) if (RunGoldStein1()) okGoldStein1++;
  StatusPrint(5, "GoldStein1", okGoldStein1, N);
  for (i = 0; i < N; i++) if (RunGoldStein2()) okGoldStein2++;
  StatusPrint(6, "GoldStein2", okGoldStein2, N);
  if (RunTrigoFletcher()) okTrigoFletcher++;
  StatusPrint(7, "TrigoFletcher", okTrigoFletcher, 1);

  gBenchmark->Stop("stressFit");
  
#if VERBOSE 
  gBenchmark->Print("stressFit");
  printf("******************************************************************\n");
#endif

#if VERBOSE
  Double_t reftime = 12.07; //macbrun compiled

  Double_t rootmarks = 800.*reftime/gBenchmark->GetCpuTime("stressFit");
  
  printf("******************************************************************\n");
  printf("* stressFit * ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
         gROOT->GetVersionDate(),gROOT->GetVersionTime());
  printf("******************************************************************\n");
#endif

  return gBenchmark->GetCpuTime("stressFit");
}

Double_t stressLinear(Int_t maxSizeReq,Int_t verbose)
{
#if VERBOSE
  cout << "******************************************************************" <<endl;
  cout << "*  Starting  Linear Algebra - S T R E S S suite                  *" <<endl;
  cout << "******************************************************************" <<endl;
  cout << "******************************************************************" <<endl;

  gVerbose = verbose;
#endif

  gBenchmark->Start("stressLinear");

  gNrLoop = nrSize-1;
  while (gNrLoop > 0 && maxSizeReq < gSizeA[gNrLoop])
    gNrLoop--;

  const Int_t maxSize = gSizeA[gNrLoop];

  // Matrix
  {
#if VERBOSE
    cout << "*  Starting  Matrix - S T R E S S                                *" <<endl;
    cout << "******************************************************************" <<endl;
#endif
    mstress_allocation(maxSize);
    mstress_matrix_fill(maxSize,maxSize/2);
    mstress_element_op(maxSize,maxSize/2);
    mstress_binary_ebe_op(maxSize/2,maxSize);
    mstress_transposition(maxSize);
    mstress_special_creation(maxSize);
    mstress_matrix_promises(maxSize);
    mstress_norms(maxSize,maxSize/2);
    mstress_determinant(maxSize);
    mstress_mm_multiplications();
    mstress_sym_mm_multiplications(maxSize);
    mstress_vm_multiplications();
    mstress_inversion();

    mstress_matrix_io();
#if VERBOSE    
    cout << "******************************************************************" <<endl;
#endif    
  }

  // Sparse Matrix
  {
#if VERBOSE        
    cout << "*  Starting  Sparse Matrix - S T R E S S                         *" <<endl;
    cout << "******************************************************************" <<endl;
#endif    
    spstress_allocation(maxSize);
    spstress_matrix_fill(maxSize,maxSize/2);
    spstress_element_op(maxSize,maxSize/2);
    spstress_binary_ebe_op(maxSize/2,maxSize);
    spstress_transposition(maxSize);
    spstress_norms(maxSize,maxSize/2);
    spstress_mm_multiplications();
    spstress_vm_multiplications();
    spstress_matrix_slices(maxSize);
    spstress_matrix_io();
#if VERBOSE        
    cout << "******************************************************************" <<endl;
#endif
  }

  {
#if VERBOSE    
    cout << "*  Starting  Vector - S T R E S S                                *" <<endl;
    cout << "******************************************************************" <<endl;
#endif    
    vstress_allocation(maxSize);
    vstress_element_op(maxSize);
    vstress_binary_op(maxSize);
    vstress_norms(maxSize);
    vstress_matrix_slices(maxSize);
    vstress_vector_io();
#if VERBOSE        
    cout << "******************************************************************" <<endl;
#endif
  }

  // Linear Algebra
  {
#if VERBOSE    
    cout << "*  Starting  Linear Algebra - S T R E S S                        *" <<endl;
    cout << "******************************************************************" <<endl;
#endif
    astress_decomp();
    astress_lineqn();
    astress_pseudo();
    astress_eigen(5);
    astress_decomp_io(10);
#if VERBOSE    
    cout << "******************************************************************" <<endl;
#endif    
  }

  //Backward Compatibility of Streamers
  {
#if VERBOSE
    cout << "*  Starting  Backward IO compatibility - S T R E S S             *" <<endl;
    cout << "******************************************************************" <<endl;
#endif    
    stress_backward_io();
#if VERBOSE
    cout << "******************************************************************" <<endl;
#endif
  }

  gBenchmark->Stop("stressLinear");

#if VERBOSE
  printf("******************************************************************\n");
  gBenchmark->Print("stressLinear");
#endif

  const Int_t nr = 7;
  const Double_t x_b12[] = { 10.,   30.,   50.,   100.,  300.,  500.,    700.};
  const Double_t y_b12[] = {10.74, 15.72, 20.00, 35.79, 98.77, 415.34, 1390.33};

  TGraph gr(nr,x_b12,y_b12);
  Double_t ct = gBenchmark->GetCpuTime("stressLinear");

#if VERBOSE
  Double_t rootmarks = 600*gr.Eval(maxSize)/ct;
  printf("******************************************************************\n");
  printf("*  stressLinear * ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
         gROOT->GetVersionDate(),gROOT->GetVersionTime());
  printf("******************************************************************\n");
#endif

  return ct;
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void mstress_allocation(Int_t msize)
{
  if (gVerbose)
    cout << "\n\n---> Test allocation and compatibility check" << endl;

  Int_t i,j;
  Bool_t ok = kTRUE;

  TMatrixD m1(4,msize);
  for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
    for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
      m1(i,j) = TMath::Pi()*i+TMath::E()*j;

  TMatrixD m2(0,3,0,msize-1);
  TMatrixD m3(1,4,0,msize-1);
  TMatrixD m4(m1);

  if (gVerbose) {
    cout << "\nStatus information reported for matrix m3:" << endl;
    cout << "  Row lower bound ... " << m3.GetRowLwb() << endl;
    cout << "  Row upper bound ... " << m3.GetRowUpb() << endl;
    cout << "  Col lower bound ... " << m3.GetColLwb() << endl;
    cout << "  Col upper bound ... " << m3.GetColUpb() << endl;
    cout << "  No. rows ..........." << m3.GetNrows()  << endl;
    cout << "  No. cols ..........." << m3.GetNcols()  << endl;
    cout << "  No. of elements ...." << m3.GetNoElements() << endl;
  }

  if (gVerbose)
    cout << "\nCheck matrices 1 & 2 for compatibility" << endl;
  ok &= AreCompatible(m1,m2,gVerbose);

  if (gVerbose)
    cout << "Check matrices 1 & 4 for compatibility" << endl;
  ok &= AreCompatible(m1,m4,gVerbose);

  if (gVerbose)
    cout << "m2 has to be compatible with m3 after resizing to m3" << endl;
  m2.ResizeTo(m3);
  ok &= AreCompatible(m2,m3,gVerbose);

  TMatrixD m5(m1.GetNrows()+1,m1.GetNcols()+5);
  for (i = m5.GetRowLwb(); i <= m5.GetRowUpb(); i++)
    for (j = m5.GetColLwb(); j <= m5.GetColUpb(); j++)
      m5(i,j) = TMath::Pi()*i+TMath::E()*j;

  if (gVerbose)
    cout << "m1 has to be compatible with m5 after resizing to m5" << endl;
  m1.ResizeTo(m5.GetNrows(),m5.GetNcols());
  ok &= AreCompatible(m1,m5,gVerbose);

  if (gVerbose)
    cout << "m1 has to be equal to m4 after stretching and shrinking" << endl;
  m1.ResizeTo(m4.GetNrows(),m4.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m4,gVerbose,EPSILON);
  if (gVerbose)
    cout << "m5 has to be equal to m1 after shrinking" << endl;
  m5.ResizeTo(m1.GetNrows(),m1.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m5,gVerbose,EPSILON);

  if (gVerbose)
    cout << "stretching and shrinking for small matrices (stack)" << endl;
  if (gVerbose)
    cout << "m8 has to be equal to m7 after stretching and shrinking" << endl;
  TMatrixD m6(4,4);
  for (i = m6.GetRowLwb(); i <= m6.GetRowUpb(); i++)
    for (j = m6.GetColLwb(); j <= m6.GetColUpb(); j++)
      m6(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixD m8(3,3);
  for (i = m8.GetRowLwb(); i <= m8.GetRowUpb(); i++)
    for (j = m8.GetColLwb(); j <= m8.GetColUpb(); j++)
      m8(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixD m7(m8);

  m8.ResizeTo(4,4);
  m8.ResizeTo(3,3);
  ok &= VerifyMatrixIdentity(m7,m8,gVerbose,EPSILON);

  if (gVerbose)
    cout << "m6 has to be equal to m8 after shrinking" << endl;
  m6.ResizeTo(3,3);
  ok &= VerifyMatrixIdentity(m6,m8,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(1,"Allocation, Resizing",ok);
}

class FillMatrix : public TElementPosActionD {
   Int_t no_elems,no_cols;
   void Operation(Double_t &element) const
      { element = 4*TMath::Pi()/no_elems * (fI*no_cols+fJ); }
public:
   FillMatrix() {}
   FillMatrix(const TMatrixD &m) :
         no_elems(m.GetNoElements()),no_cols(m.GetNcols()) { }
};

//
//------------------------------------------------------------------------
//          Test Filling of matrix
//
void mstress_matrix_fill(Int_t rsize,Int_t csize)
{
  if (gVerbose)
    cout << "\n\n---> Test different matrix filling methods\n" << endl;

  Bool_t ok = kTRUE;
  if (gVerbose)
    cout << "Creating m  with Apply function..." << endl;
  TMatrixD m(-1,rsize-2,1,csize);
#ifndef __CINT__
  FillMatrix f(m);
  m.Apply(f);
#else
  for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
    for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
      m(i,j) = 4*TMath::Pi()/m.GetNoElements() * (i*m.GetNcols()+j);
#endif

  {
    if (gVerbose)
      cout << "Check identity between m and matrix filled through (i,j)" << endl;

    TMatrixD m_overload1(-1,rsize-2,1,csize);
    TMatrixD m_overload2(-1,rsize-2,1,csize);

    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
    {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
      {
        const Double_t val = 4*TMath::Pi()/rsize/csize*(i*csize+j);
        m_overload1(i,j)  = val;
        m_overload2[i][j] = val;
      }
    }

    ok &= VerifyMatrixIdentity(m,m_overload1,gVerbose,EPSILON);
    if (gVerbose)
      cout << "Check identity between m and matrix filled through [i][j]" << endl;
    ok &= VerifyMatrixIdentity(m,m_overload2,gVerbose,EPSILON);
    if (gVerbose)
      cout << "Check identity between matrix filled through [i][j] and (i,j)" << endl;
    ok &= VerifyMatrixIdentity(m_overload1,m_overload2,gVerbose,EPSILON);
  }

  {
    TArrayD a_fortran(rsize*csize);
    TArrayD a_c      (rsize*csize);
    for (Int_t i = 0; i < rsize; i++)
    {
      for (Int_t j = 0; j < csize; j++)
      {
        a_c[i*csize+j]       = 4*TMath::Pi()/rsize/csize*((i-1)*csize+j+1);
        a_fortran[i+rsize*j] = a_c[i*csize+j];
      }
    }

    if (gVerbose)
      cout << "Creating m_fortran by filling with fortran stored matrix" << endl;
    TMatrixD m_fortran(-1,rsize-2,1,csize,a_fortran.GetArray(),"F");
    if (gVerbose)
      cout << "Check identity between m and m_fortran" << endl;
    ok &= VerifyMatrixIdentity(m,m_fortran,gVerbose,EPSILON);

    if (gVerbose)
      cout << "Creating m_c by filling with c stored matrix" << endl;
    TMatrixD m_c(-1,rsize-2,1,csize,a_c.GetArray());
    if (gVerbose)
      cout << "Check identity between m and m_c" << endl;
    ok &= VerifyMatrixIdentity(m,m_c,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      cout << "Check insertion/extraction of sub-matrices" << endl;
    {
      TMatrixD m_sub1 = m;
      m_sub1.ResizeTo(0,rsize-2,2,csize);
      TMatrixD m_sub2 = m.GetSub(0,rsize-2,2,csize,"");
      ok &= VerifyMatrixIdentity(m_sub1,m_sub2,gVerbose,EPSILON);
    }

    {
      TMatrixD m2(-1,rsize-2,1,csize);
      TMatrixD m_part1 = m.GetSub(0,rsize-2,2,csize,"");
      TMatrixD m_part2 = m.GetSub(0,rsize-2,1,1,"");
      TMatrixD m_part3 = m.GetSub(-1,-1,2,csize,"");
      TMatrixD m_part4 = m.GetSub(-1,-1,1,1,"");
      m2.SetSub(0,2,m_part1);
      m2.SetSub(0,1,m_part2);
      m2.SetSub(-1,2,m_part3);
      m2.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m2,gVerbose,EPSILON);
    }

    {
      TMatrixD m2(-1,rsize-2,1,csize);
      TMatrixD m_part1 = m.GetSub(0,rsize-2,2,csize,"S");
      TMatrixD m_part2 = m.GetSub(0,rsize-2,1,1,"S");
      TMatrixD m_part3 = m.GetSub(-1,-1,2,csize,"S");
      TMatrixD m_part4 = m.GetSub(-1,-1,1,1,"S");
      m2.SetSub(0,2,m_part1);
      m2.SetSub(0,1,m_part2);
      m2.SetSub(-1,2,m_part3);
      m2.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m2,gVerbose,EPSILON);
    }
  }

  {
    if (gVerbose)
      cout << "Check sub-matrix views" << endl;
    {
      TMatrixD m3(-1,rsize-2,1,csize);
      TMatrixDSub(m3,0,rsize-2,2,csize) = TMatrixDSub(m,0,rsize-2,2,csize);
      TMatrixDSub(m3,0,rsize-2,1,1)     = TMatrixDSub(m,0,rsize-2,1,1);
      TMatrixDSub(m3,-1,-1,2,csize)     = TMatrixDSub(m,-1,-1,2,csize);
      TMatrixDSub(m3,-1,-1,1,1)         = TMatrixDSub(m,-1,-1,1,1);
      ok &= VerifyMatrixIdentity(m,m3,gVerbose,EPSILON);

      TMatrixD unit(3,3);
      TMatrixDSub(m3,1,3,1,3)  = unit.UnitMatrix();
      TMatrixDSub(m3,1,3,1,3) *= m.GetSub(1,3,1,3);
      ok &= VerifyMatrixIdentity(m,m3,gVerbose,EPSILON);

      TMatrixDSub(m3,0,rsize-2,2,csize) = 1.0;
      TMatrixDSub(m3,0,rsize-2,1,1)     = 1.0;
      TMatrixDSub(m3,-1,-1,2,csize)     = 1.0;
      TMatrixDSub(m3,-1,-1,1,1)         = 1.0;
      ok &= (m3 == 1.0);

    }
  }

  {
    if (gVerbose)
      cout << "Check array Use" << endl;
    {
      TMatrixD *m1a = new TMatrixD(m);
      TMatrixD *m2a = new TMatrixD();
      m2a->Use(m1a->GetRowLwb(),m1a->GetRowUpb(),m1a->GetColLwb(),m1a->GetColUpb(),m1a->GetMatrixArray());
      ok &= VerifyMatrixIdentity(m,*m2a,gVerbose,EPSILON);
      m2a->Sqr();
      TMatrixD m4 = m; m4.Sqr();
      ok &= VerifyMatrixIdentity(m4,*m1a,gVerbose,EPSILON);
      delete m1a;
      delete m2a;
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(2,"Filling, Inserting, Using",ok);
}

//
//------------------------------------------------------------------------
//                Test uniform element operations
//
typedef  Double_t (*dfunc)(Double_t);
class ApplyFunction : public TElementActionD {
   dfunc fFunc;
   void Operation(Double_t &element) const { element = fFunc(Double_t(element)); }
public:
   ApplyFunction(dfunc func) : fFunc(func) { }
};

void mstress_element_op(Int_t rsize,Int_t csize)
{
  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TMatrixD m(-1,rsize-2,1,csize);

  if (gVerbose)
    cout << "\nWriting zeros to m..." << endl;
  {
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for(Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = 0;
    ok &= VerifyMatrixValue(m,0.,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "Creating zero m1 ..." << endl;
  TMatrixD m1(TMatrixD::kZero, m);
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Comparing m1 with 0 ..." << endl;
  R__ASSERT(m1 == 0);
  R__ASSERT(!(m1 != 0));

  if (gVerbose)
    cout << "Writing a pattern " << pattern << " by assigning to m(i,j)..." << endl;
  {
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = pattern;
    ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "Writing the pattern by assigning to m1 as a whole ..."  << endl;
  m1 = pattern;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Comparing m and m1 ..." << endl;
  R__ASSERT(m == m1);
  if (gVerbose)
    cout << "Comparing (m=0) and m1 ..." << endl;
  R__ASSERT(!(m.Zero() == m1));

  if (gVerbose)
    cout << "Clearing m1 ..." << endl;
  m1.Zero();
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nClear m and add the pattern" << endl;
  m.Zero();
  m += pattern;
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   add the doubled pattern with the negative sign" << endl;
  m += -2*pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   subtract the trippled pattern with the negative sign" << endl;
  m -= -3*pattern;
  ok &= VerifyMatrixValue(m,2*pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nVerify comparison operations when all elems are the same" << endl;
  m = pattern;
  R__ASSERT( m == pattern && !(m != pattern) );
  R__ASSERT( m > 0 && m >= pattern && m <= pattern );
  R__ASSERT( m > -pattern && m >= -pattern );
  R__ASSERT( m <= pattern && !(m < pattern) );
  m -= 2*pattern;
  R__ASSERT( m  < -pattern/2 && m <= -pattern/2 );
  R__ASSERT( m  >= -pattern && !(m > -pattern) );

  if (gVerbose)
    cout << "\nVerify comparison operations when not all elems are the same" << endl;
  m = pattern; m(m.GetRowUpb(),m.GetColUpb()) = pattern-1;
  R__ASSERT( !(m == pattern) && !(m != pattern) );
  R__ASSERT( m != 0 );                   // none of elements are 0
  R__ASSERT( !(m >= pattern) && m <= pattern && !(m<pattern) );
  R__ASSERT( !(m <= pattern-1) && m >= pattern-1 && !(m>pattern-1) );

  if (gVerbose)
    cout << "\nAssign 2*pattern to m by repeating additions" << endl;
  m = 0; m += pattern; m += pattern;
  if (gVerbose)
    cout << "Assign 2*pattern to m1 by multiplying by two " << endl;
  m1 = pattern; m1 *= 2;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    cout << "Multiply m1 by one half returning it to the 1*pattern" << endl;
  m1 *= 1/2.;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nAssign -pattern to m and m1" << endl;
  m.Zero(); m -= pattern; m1 = -pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    cout << "m = sqrt(sqr(m)); m1 = abs(m1); Now m and m1 have to be the same" << endl;
  m.Sqr();
  ok &= VerifyMatrixValue(m,pattern*pattern,gVerbose,EPSILON);
  m.Sqrt();
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  m1.Abs();
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  ok &= VerifyMatrixIdentity(m1,m,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nCheck out to see that sin^2(x) + cos^2(x) = 1" << endl;
  {
#ifndef __CINT__
    FillMatrix f(m);
    m.Apply(f);
#else
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = 4*TMath::Pi()/m.GetNoElements() * (i*m.GetNcols()+j);
#endif
  }
  m1 = m;
  {
#ifndef __CINT__
    ApplyFunction s(&TMath::Sin);
    ApplyFunction c(&TMath::Cos);
    m.Apply(s);
    m1.Apply(c);
#else
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        m(i,j)  = TMath::Sin(m(i,j));
        m1(i,j) = TMath::Cos(m1(i,j));
      }
    }
#endif
  }
  m.Sqr();
  m1.Sqr();
  m += m1;
  ok &= VerifyMatrixValue(m,1.,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(3,"Uniform matrix operations",ok);
}

//
//------------------------------------------------------------------------
//        Test binary matrix element-by-element operations
//
void mstress_binary_ebe_op(Int_t rsize, Int_t csize)
{
  if (gVerbose)
    cout << "\n---> Test Binary Matrix element-by-element operations" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 4.25;

  TMatrixD m(2,rsize+1,0,csize-1);
  TMatrixD m1(TMatrixD::kZero,m);
  TMatrixD mp(TMatrixD::kZero,m);

  {
    for (Int_t i = mp.GetRowLwb(); i <= mp.GetRowUpb(); i++)
      for (Int_t j = mp.GetColLwb(); j <= mp.GetColUpb(); j++)
        mp(i,j) = (i-m.GetNrows()/2.)*j*pattern;
  }

  if (gVerbose)
    cout << "\nVerify assignment of a matrix to the matrix" << endl;
  m = pattern;
  m1.Zero();
  m1 = m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  R__ASSERT( m1 == m );

  if (gVerbose)
    cout << "\nAdding the matrix to itself, uniform pattern " << pattern << endl;
  m.Zero(); m = pattern;
  m1 = m; m1 += m1;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  subtracting two matrices ..." << endl;
  m1 -= m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  subtracting the matrix from itself" << endl;
  m1 -= m1;
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  adding two matrices together" << endl;
  m1 += m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose) {
    cout << "\nArithmetic operations on matrices with not the same elements" << endl;
    cout << "   adding mp to the zero matrix..." << endl;
  }
  m.Zero(); m += mp;
  ok &= VerifyMatrixIdentity(m,mp,gVerbose,EPSILON);
  m1 = m;
  if (gVerbose)
    cout << "   making m = 3*mp and m1 = 3*mp, via add() and succesive mult" << endl;
  Add(m,2.,mp);
  m1 += m1; m1 += mp;
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   clear both m and m1, by subtracting from itself and via add()" << endl;
  m1 -= m1;
  Add(m,-3.,mp);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

  if (gVerbose) {
    cout << "\nTesting element-by-element multiplications and divisions" << endl;
    cout << "   squaring each element with sqr() and via multiplication" << endl;
  }
  m = mp; m1 = mp;
  m.Sqr();
  ElementMult(m1,m1);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   compare (m = pattern^2)/pattern with pattern" << endl;
  m = pattern; m1 = pattern;
  m.Sqr();
  ElementDiv(m,m1);
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    Compare(m1,m);

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(4,"Binary Matrix element-by-element operations",ok);
}

//
//------------------------------------------------------------------------
//              Verify matrix transposition
//
void mstress_transposition(Int_t msize)
{
  if (gVerbose) {
    cout << "\n---> Verify matrix transpose "
            "for matrices of a characteristic size " << msize << endl;
  }

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      cout << "\nCheck to see that a square UnitMatrix stays the same";
    TMatrixD m(msize,msize);
    m.UnitMatrix();
    TMatrixD mt(TMatrixD::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "\nTest a non-square UnitMatrix";
    TMatrixD m(msize,msize+1);
    m.UnitMatrix();
    TMatrixD mt(TMatrixD::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows() );
    for (Int_t i = m.GetRowLwb(); i <= TMath::Min(m.GetRowUpb(),m.GetColUpb()); i++)
      for (Int_t j = m.GetColLwb(); j <= TMath::Min(m.GetRowUpb(),m.GetColUpb()); j++)
        ok &= ( m(i,j) == mt(i,j) ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "\nCheck to see that a symmetric (Hilbert)Matrix stays the same";
    TMatrixD m = THilbertMatrixD(msize,msize);
    TMatrixD mt(TMatrixD::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "\nCheck transposing a non-symmetric matrix";
    TMatrixD m = THilbertMatrixD(msize+1,msize);
    m(1,2) = TMath::Pi();
    TMatrixD mt(TMatrixD::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows());
    R__ASSERT(mt(2,1)  == (Double_t)TMath::Pi() && mt(1,2)  != (Double_t)TMath::Pi());
    R__ASSERT(mt[2][1] == (Double_t)TMath::Pi() && mt[1][2] != (Double_t)TMath::Pi());

    if (gVerbose)
      cout << "\nCheck double transposing a non-symmetric matrix" << endl;
    TMatrixD mtt(TMatrixD::kTransposed,mt);
    ok &= ( m == mtt ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(5,"Matrix transposition",ok);
}

//
//------------------------------------------------------------------------
//           Test special matrix creation
//
class MakeHilbert : public TElementPosActionD {
  void Operation(Double_t &element) const { element = 1./(fI+fJ+1); }
public:
  MakeHilbert() { }
};

#if !defined (__CINT__) || defined (__MAKECINT__)
class TestUnit : public TElementPosActionD {
  mutable Int_t fIsUnit;
  void Operation(Double_t &element) const
      { if (fIsUnit)
          fIsUnit = ((fI==fJ) ? (element == 1.0) : (element == 0)); }
public:
  TestUnit() : fIsUnit(0==0) { }
  Int_t is_indeed_unit() const { return fIsUnit; }
};
#else
  Bool_t is_indeed_unit(TMatrixD &m) {
    Bool_t isUnit = kTRUE;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        if (isUnit)
          isUnit = ((i==j) ? (m(i,j) == 1.0) : (m(i,j) == 0));
      }
    return isUnit;
  }
#endif

void mstress_special_creation(Int_t dim)
{
  if (gVerbose)
    cout << "\n---> Check creating some special matrices of dimension " << dim << endl;

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      cout << "\ntest creating Hilbert matrices" << endl;
    TMatrixD m = THilbertMatrixD(dim+1,dim);
    TMatrixD m1(TMatrixD::kZero,m);
    ok &= ( !(m == m1) ) ? kTRUE : kFALSE;
    ok &= ( m != 0 ) ? kTRUE : kFALSE;
#ifndef __CINT__
    MakeHilbert mh;
    m1.Apply(mh);
#else
    for (Int_t i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
      for (Int_t j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
        m1(i,j) = 1./(i+j+1);
#endif
    ok &= ( m1 != 0 ) ? kTRUE : kFALSE;
    ok &= ( m == m1 ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "test creating zero matrix and copy constructor" << endl;
    TMatrixD m = THilbertMatrixD(dim,dim+1);
    ok &= ( m != 0 ) ? kTRUE : kFALSE;
    TMatrixD m1(m);               // Applying the copy constructor
    ok &= ( m1 == m ) ? kTRUE : kFALSE;
    TMatrixD m2(TMatrixD::kZero,m);
    ok &= ( m2 == 0 ) ? kTRUE : kFALSE;
    ok &= ( m != 0 ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "test creating unit matrices" << endl;
    TMatrixD m(dim,dim);
#ifndef __CINT__
    {
      TestUnit test_unit;
      m.Apply(test_unit);
      ok &= ( !test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( !is_indeed_unit(m) ) ? kTRUE : kFALSE;
#endif
    m.UnitMatrix();
#ifndef __CINT__
    {
      TestUnit test_unit;
       m.Apply(test_unit);
       ok &= ( test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( is_indeed_unit(m) ) ? kTRUE : kFALSE;
#endif
    m.ResizeTo(dim-1,dim);
    TMatrixD m2(TMatrixD::kUnit,m);
#ifndef __CINT__
    {
      TestUnit test_unit;
      m2.Apply(test_unit);
      ok &= ( test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( is_indeed_unit(m2) ) ? kTRUE : kFALSE;
#endif
    m.ResizeTo(dim,dim-2);
    m.UnitMatrix();
#ifndef __CINT__
    {
      TestUnit test_unit;
      m.Apply(test_unit);
      ok &= ( test_unit.is_indeed_unit() ) ? kTRUE : kFALSE;
    }
#else
    ok &= ( is_indeed_unit(m) ) ? kTRUE : kFALSE;
#endif
  }

  {
    if (gVerbose)
      cout << "check to see that Haar matrix has *exactly* orthogonal columns" << endl;
    Int_t j;
    const Int_t order = 5;
    const TMatrixD haar = THaarMatrixD(order);
    ok &= ( haar.GetNrows() == (1<<order) &&
               haar.GetNrows() == haar.GetNcols() ) ? kTRUE : kFALSE;
    TVectorD colj(1<<order);
    TVectorD coll(1<<order);
    for (j = haar.GetColLwb(); j <= haar.GetColUpb(); j++) {
      colj = TMatrixDColumn_const(haar,j);
      ok &= (TMath::Abs(colj*colj-1.0) <= 1.0e-15 ) ? kTRUE : kFALSE;
      for (Int_t l = j+1; l <= haar.GetColUpb(); l++) {
        coll = TMatrixDColumn_const(haar,l);
        const Double_t val = colj*coll;
        ok &= ( TMath::Abs(val) <= 1.0e-15 ) ? kTRUE : kFALSE;
      }
    }

    if (gVerbose)
      cout << "make Haar (sub)matrix and test it *is* a submatrix" << endl;
    const Int_t no_sub_cols = (1<<order) - 3;
    const TMatrixD haar_sub = THaarMatrixD(order,no_sub_cols);
    ok &= ( haar_sub.GetNrows() == (1<<order) &&
               haar_sub.GetNcols() == no_sub_cols ) ? kTRUE : kFALSE;
    for (j = haar_sub.GetColLwb(); j <= haar_sub.GetColUpb(); j++) {
      colj = TMatrixDColumn_const(haar,j);
      coll = TMatrixDColumn_const(haar_sub,j);
      ok &= VerifyVectorIdentity(colj,coll,gVerbose);
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(6,"Haar/Hilbert Matrix",ok);
}

#ifndef __CINT__
//
//------------------------------------------------------------------------
//           Test matrix promises
//
class hilbert_matrix_promise : public TMatrixDLazy {
  void FillIn(TMatrixD &m) const { m = THilbertMatrixD(m.GetRowLwb(),m.GetRowUpb(),
                                                   m.GetColLwb(),m.GetColUpb()); }

public:
  hilbert_matrix_promise(Int_t nrows,Int_t ncols)
     : TMatrixDLazy(nrows,ncols) {}
  hilbert_matrix_promise(Int_t row_lwb,Int_t row_upb,
                         Int_t col_lwb,Int_t col_upb)
     : TMatrixDLazy(row_lwb,row_upb,col_lwb,col_upb) { }
};

void mstress_matrix_promises(Int_t dim)
{
  if (gVerbose)
    cout << "\n---> Check making/forcing promises, (lazy)matrices of dimension " << dim << endl;

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      cout << "\nmake a promise and force it by a constructor" << endl;
    TMatrixD m  = hilbert_matrix_promise(dim,dim+1);
    TMatrixD m1 = THilbertMatrixD(dim,dim+1);
    ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      cout << "make a promise and force it by an assignment" << endl;
    TMatrixD m(-1,dim,0,dim);
    m = hilbert_matrix_promise(-1,dim,0,dim);
    TMatrixD m1 = THilbertMatrixD(-1,dim,0,dim);
    ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(7,"Matrix promises",ok);
}
#endif

//
//------------------------------------------------------------------------
//             Verify the norm calculation
//
void mstress_norms(Int_t rsize_req,Int_t csize_req)
{
  if (gVerbose)
    cout << "\n---> Verify norm calculations" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 10.25;

  Int_t rsize = rsize_req;
  if (rsize%2 != 0)  rsize--;
  Int_t csize = csize_req;
  if (csize%2 != 0)  csize--;
  if (rsize%2 == 1 || csize%2 == 1) {
    cout << "rsize: " << rsize <<endl;
    cout << "csize: " << csize <<endl;
    Fatal("mstress_norms","Sorry, size of the matrix to test must be even for this test\n");
  }

  TMatrixD m(rsize,csize);

  if (gVerbose)
    cout << "\nAssign " << pattern << " to all the elements and check norms" << endl;
  m = pattern;
  if (gVerbose)
    cout << "  1. (col) norm should be pattern*nrows" << endl;
  ok &= ( m.Norm1() == pattern*m.GetNrows() ) ? kTRUE : kFALSE;
  ok &= ( m.Norm1() == m.ColNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Inf (row) norm should be pattern*ncols" << endl;
  ok &= ( m.NormInf() == pattern*m.GetNcols() ) ? kTRUE : kFALSE;
  ok &= ( m.NormInf() == m.RowNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Square of the Eucl norm has got to be pattern^2 * no_elems" << endl;
  ok &= ( m.E2Norm() == (pattern*pattern)*m.GetNoElements() ) ? kTRUE : kFALSE;
  TMatrixD m1(TMatrixD::kZero,m);
  ok &= ( m.E2Norm() == E2Norm(m,m1) ) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(8,"Matrix Norms",ok);
}

//
//------------------------------------------------------------------------
//              Verify the determinant evaluation
//
void mstress_determinant(Int_t msize)
{
  if (gVerbose)
    cout << "\n---> Verify determinant evaluation for a square matrix of size " << msize << endl;

  Bool_t ok = kTRUE;
  TMatrixD m(msize,msize);
  const Double_t pattern = 2.5;

  {
    if (gVerbose)
      cout << "\nCheck to see that the determinant of the unit matrix is one" <<endl;
    m.UnitMatrix();
    Double_t d1,d2;
    m.Determinant(d1,d2);
    const Double_t det = d1*TMath::Power(2.,d2);
    if (gVerbose) {
      cout << "det = " << det << " deviation= " << TMath::Abs(det-1);
      cout << ( (TMath::Abs(det-1) < EPSILON) ? " OK" : " too large") <<endl;
    }
    ok &= (TMath::Abs(det-1) < EPSILON) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      cout << "\nCheck the determinant for the matrix with " << pattern << " at the diagonal" << endl;
    {
      for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
        for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
          m(i,j) = ( i==j ? pattern : 0 );
    }
    Double_t d1,d2;
    m.Determinant(d1,d2);
    const Double_t d1_abs = TMath::Abs(d1);
    const Double_t log_det1 = d2*TMath::Log(2.)+TMath::Log(d1_abs);
    const Double_t log_det2 = ((Double_t)m.GetNrows())*TMath::Log(pattern);
    if (gVerbose) {
      cout << "log_det1 = " << log_det1 << " deviation= " << TMath::Abs(log_det1-log_det2);
      cout << ( (TMath::Abs(log_det1-log_det2) < msize*EPSILON) ? " OK" : " too large") <<endl;
    }
    ok &= ( TMath::Abs(log_det1-log_det2) < msize*EPSILON  ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      cout << "\nCheck the determinant of the transposed matrix" << endl;
    m.UnitMatrix();
    m(1,2) = pattern;
    TMatrixD m_tran(TMatrixD::kTransposed,m);
    ok &= ( !(m == m_tran) ) ? kTRUE : kFALSE;
    const Double_t det1 = m.Determinant();
    const Double_t det2 = m_tran.Determinant();
    if (gVerbose) {
      cout << "det1 = " << det1 << " deviation= " << TMath::Abs(det1-det2);
      cout << ( (TMath::Abs(det1-det2) < msize*EPSILON) ? " OK" : " too large") <<endl;
    }
    ok &= ( TMath::Abs(det1-det2) < msize*EPSILON ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      cout << "\nswap two rows/cols of a matrix through method 1 and watch det's sign" << endl;
    m.UnitMatrix();
    TMatrixDRow(m,3) = pattern;
    Double_t d1,d2;
    m.Determinant(d1,d2);
    TMatrixDRow row1(m,1);
    TVectorD vrow1(m.GetRowLwb(),m.GetRowUpb()); vrow1 = row1;
    TVectorD vrow3(m.GetRowLwb(),m.GetRowUpb()); vrow3 = TMatrixDRow(m,3);
    row1 = vrow3; TMatrixDRow(m,3) = vrow1;
    Double_t d1_p,d2_p;
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1+d1_p);
      cout << ( ( d1 == -d1_p ) ? " OK" : " too large") <<endl;
    }
    ok &= ( d1 == -d1_p ) ? kTRUE : kFALSE;
    TMatrixDColumn col2(m,2);
    TVectorD vcol2(m.GetRowLwb(),m.GetRowUpb()); vcol2 = col2;
    TVectorD vcol4(m.GetRowLwb(),m.GetRowUpb()); vcol4 = TMatrixDColumn(m,4);
    col2 = vcol4; TMatrixDColumn(m,4) = vcol2;
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1-d1_p);
      cout << ( ( d1 == d1_p ) ? " OK" : " too large") <<endl;
    }
    ok &= ( d1 == d1_p ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      cout << "\nswap two rows/cols of a matrix through method 2 and watch det's sign" << endl;
    m.UnitMatrix();
    TMatrixDRow(m,3) = pattern;
    Double_t d1,d2;
    m.Determinant(d1,d2);

    TMatrixD m_save( m);
    TMatrixDRow(m,1) = TMatrixDRow(m_save,3);
    TMatrixDRow(m,3) = TMatrixDRow(m_save,1);
    Double_t d1_p,d2_p;
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1+d1_p);
      cout << ( ( d1 == -d1_p ) ? " OK" : " too large") <<endl;
    }
    ok &= ( d1 == -d1_p ) ? kTRUE : kFALSE;

    m_save = m;
    TMatrixDColumn(m,2) = TMatrixDColumn(m_save,4);
    TMatrixDColumn(m,4) = TMatrixDColumn(m_save,2);
    m.Determinant(d1_p,d2_p);
    if (gVerbose) {
      cout << "d1 = " << d1 << " deviation= " << TMath::Abs(d1-d1_p);
      cout << ( ( d1 == d1_p ) ? " OK" : " too large") <<endl;
    }
    ok &= ( d1 == d1_p ) ? kTRUE : kFALSE;
  }

  if (ok)
  {
    if (gVerbose)
      cout << "\nCheck the determinant for the matrix with " << pattern << " at the anti-diagonal" << endl;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++)
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        m(i,j) = ( i==(m.GetColUpb()+m.GetColLwb()-j) ? pattern : 0 );

    Double_t d1,d2;
    m.Determinant(d1,d2);
    const Double_t d1_abs = TMath::Abs(d1);
    const Double_t log_det1 = d2*TMath::Log(2.)+TMath::Log(d1_abs);
    const Double_t log_det2 = ((Double_t)m.GetNrows())*TMath::Log(pattern);
    const Double_t sign     = ( m.GetNrows()*(m.GetNrows()-1)/2 & 1 ? -1 : 1 );
    if (gVerbose) {
      cout << "log_det1 = " << log_det1 << " deviation= " << TMath::Abs(log_det1-log_det2);
      cout << ( (TMath::Abs(log_det1-log_det2) < msize*EPSILON) ? " OK" : " too large") <<endl;
      cout << ( sign * d1 > 0. ? " OK sign" : " wrong sign") <<endl;
    }
    ok &= ( TMath::Abs(log_det1-log_det2) < msize*EPSILON  ) ? kTRUE : kFALSE;
    ok &= ( sign * d1 > 0. ) ? kTRUE : kFALSE;
  }

  // Next test disabled because it produces (of course) a Warning
  if (0) {
    if (gVerbose)
      cout << "\nCheck the determinant for the singular matrix"
              "\n\tdefined as above with zero first row" << endl;
    m.Zero();
    {
      for (Int_t i = m.GetRowLwb()+1; i <= m.GetRowUpb(); i++)
        for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++)
          m(i,j) = ( i==(m.GetColUpb()+m.GetColLwb()-j) ? pattern : 0 );
    }
    cout << "\n\tdeterminant is " << m.Determinant();
    ok &= ( m.Determinant() == 0 ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nCheck out the determinant of the Hilbert matrix";
  TMatrixDSym H = THilbertMatrixDSym(3);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    cout << "\n    3x3 Hilbert matrix: exact determinant 1/2160 ";
    cout << "\n                              computed    1/"<< 1/H.Determinant();
  }

  H.ResizeTo(4,4);
  H = THilbertMatrixDSym(4);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    cout << "\n    4x4 Hilbert matrix: exact determinant 1/6048000 ";
    cout << "\n                              computed    1/"<< 1/H.Determinant();
  }

  H.ResizeTo(5,5);
  H = THilbertMatrixDSym(5);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    cout << "\n    5x5 Hilbert matrix: exact determinant 3.749295e-12";
    cout << "\n                              computed    "<< H.Determinant();
  }

  if (gVerbose) {
    TDecompQRH qrh(H);
    Double_t d1,d2;
    qrh.Det(d1,d2);
    cout  << "\n qrh det = " << d1*TMath::Power(2.0,d2) <<endl;
  }

  if (gVerbose) {
    TDecompChol chol(H);
    Double_t d1,d2;
    chol.Det(d1,d2);
    cout  << "\n chol det = " << d1*TMath::Power(2.0,d2) <<endl;
  }

  if (gVerbose) {
    TDecompSVD svd(H);
    Double_t d1,d2;
    svd.Det(d1,d2);
    cout  << "\n svd det = " << d1*TMath::Power(2.0,d2) <<endl;
  }

  H.ResizeTo(7,7);
  H = THilbertMatrixDSym(7);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    cout << "\n    7x7 Hilbert matrix: exact determinant 4.8358e-25";
    cout << "\n                              computed    "<< H.Determinant();
  }

  H.ResizeTo(9,9);
  H = THilbertMatrixDSym(9);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    cout << "\n    9x9 Hilbert matrix: exact determinant 9.72023e-43";
    cout << "\n                              computed    "<< H.Determinant();
  }

  H.ResizeTo(10,10);
  H = THilbertMatrixDSym(10);
  H.SetTol(1.0e-20);
  if (gVerbose) {
    cout << "\n    10x10 Hilbert matrix: exact determinant 2.16418e-53";
    cout << "\n                              computed    "<< H.Determinant();
  }

  if (gVerbose)
  cout << "\nDone\n" << endl;

  StatusPrint(9,"Matrix Determinant",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix multiplications
//
void mstress_mm_multiplications()
{
  Bool_t ok = kTRUE;

  Int_t iloop = 0;
  Int_t nr    = 0;
  while (iloop <= gNrLoop) {
    const Int_t msize = gSizeA[iloop];
    const Double_t epsilon = EPSILON*msize/100;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    Int_t i,j;
    if (msize <= 5) {
       iloop++;
       continue;  // some references to m(3,..)
    }

    if (verbose)
      cout << "\n---> Verify matrix multiplications "
              "for matrices of the characteristic size " << msize << endl;

    {
      if (verbose)
        cout << "\nTest inline multiplications of the UnitMatrix" << endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixD u(TMatrixD::kUnit,m);
      m(3,1) = TMath::Pi();
      u *= m;
      ok &= VerifyMatrixIdentity(u,m,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Test inline multiplications by a DiagMat" << endl;
      TMatrixD m = THilbertMatrixD(msize+3,msize);
      m(1,3) = TMath::Pi();
      TVectorD v(msize);
      for (i = v.GetLwb(); i <= v.GetUpb(); i++)
        v(i) = 1+i;
      TMatrixD diag(msize,msize);
      TMatrixDDiag d(diag);
      d = v;
      TMatrixD eth = m;
      for (i = eth.GetRowLwb(); i <= eth.GetRowUpb(); i++)
        for (j = eth.GetColLwb(); j <= eth.GetColUpb(); j++)
          eth(i,j) *= v(j);
      m *= diag;
      ok &= VerifyMatrixIdentity(m,eth,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Test XPP = X where P is a permutation matrix" << endl;
      TMatrixD m = THilbertMatrixD(msize-1,msize);
      m(2,3) = TMath::Pi();
      TMatrixD eth = m;
      TMatrixD p(msize,msize);
      for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++)
        p(p.GetRowUpb()+p.GetRowLwb()-i,i) = 1;
      m *= p;
      m *= p;
      ok &= VerifyMatrixIdentity(m,eth,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Test general matrix multiplication through inline mult" << endl;
      TMatrixD m = THilbertMatrixD(msize-2,msize);
      m(3,3) = TMath::Pi();
      TMatrixD mt(TMatrixD::kTransposed,m);
      TMatrixD p = THilbertMatrixD(msize,msize);
      TMatrixDDiag(p) += 1;
      TMatrixD mp(m,TMatrixD::kMult,p);
      TMatrixD m1 = m;
      m *= p;
      ok &= VerifyMatrixIdentity(m,mp,verbose,epsilon);
      TMatrixD mp1(mt,TMatrixD::kTransposeMult,p);
      VerifyMatrixIdentity(m,mp1,verbose,epsilon);
      ok &= ( !(m1 == m) );
      TMatrixD mp2(TMatrixD::kZero,m1);
      ok &= ( mp2 == 0 );
      mp2.Mult(m1,p);
      ok &= VerifyMatrixIdentity(m,mp2,verbose,epsilon);

      if (verbose)
        cout << "Test XP=X*P  vs XP=X;XP*=P" << endl;
      TMatrixD mp3 = m1*p;
      ok &= VerifyMatrixIdentity(m,mp3,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Check n * m  == n * m_sym; m_sym * n == m * n; m_sym * m_sym == m * m" <<endl;

      const TMatrixD     n     = THilbertMatrixD(0,msize-1,0,msize-1);
      const TMatrixD     m     = n;
      const TMatrixDSym  m_sym = THilbertMatrixDSym(0,msize-1);

      const TMatrixD nm1 = n * m_sym;
      const TMatrixD nm2 = n * m;
      const TMatrixD mn1 = m_sym * n;
      const TMatrixD mn2 = m * n;
      const TMatrixD mm1 = m_sym * m_sym;
      const TMatrixD mm2 = m * m;

      ok &= VerifyMatrixIdentity(nm1,nm2,verbose,epsilon);
      ok &= VerifyMatrixIdentity(mn1,mn2,verbose,epsilon);
      ok &= VerifyMatrixIdentity(mm1,mm2,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop++;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        cout << "General Matrix Multiplications failed for size " << msize << endl;
      break;
    }
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Check to see UU' = U'U = E when U is the Haar matrix" << endl;
    const Int_t order = 5;
    const Int_t no_sub_cols = (1<<order)-5;
    TMatrixD haar_sub = THaarMatrixD(5,no_sub_cols);
    TMatrixD haar_sub_t(TMatrixD::kTransposed,haar_sub);
    TMatrixD hsths(haar_sub_t,TMatrixD::kMult,haar_sub);
    TMatrixD hsths1(TMatrixD::kZero,hsths); hsths1.Mult(haar_sub_t,haar_sub);
    TMatrixD hsths_eth(TMatrixD::kUnit,hsths);
    ok &= ( hsths.GetNrows() == no_sub_cols && hsths.GetNcols() == no_sub_cols );
    ok &= VerifyMatrixIdentity(hsths,hsths_eth,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(hsths1,hsths_eth,gVerbose,EPSILON);
    TMatrixD haar = THaarMatrixD(order);
    TMatrixD unit(TMatrixD::kUnit,haar);
    TMatrixD haar_t(TMatrixD::kTransposed,haar);
    TMatrixD hth(haar,TMatrixD::kTransposeMult,haar);
    TMatrixD hht(haar,TMatrixD::kMult,haar_t);
    TMatrixD hht1 = haar; hht1 *= haar_t;
    TMatrixD hht2(TMatrixD::kZero,haar); hht2.Mult(haar,haar_t);
    ok &= VerifyMatrixIdentity(unit,hth,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(unit,hht,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(unit,hht1,gVerbose,EPSILON);
    ok &= VerifyMatrixIdentity(unit,hht2,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(10,"General Matrix Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify symmetric matrix multiplications
//
void mstress_sym_mm_multiplications(Int_t msize)
{
  if (gVerbose)
    cout << "\n---> Verify symmetric matrix multiplications "
            "for matrices of the characteristic size " << msize << endl;

  Int_t i,j;
  Bool_t ok = kTRUE;

  const Double_t epsilon = EPSILON*msize/100;

  {
    if (gVerbose)
      cout << "\nTest inline multiplications of the UnitMatrix" << endl;
    TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
    TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
    TMatrixDSym u(TMatrixDSym::kUnit,m_sym);
    TMatrixD u2 = u * m_sym;
    ok &= VerifyMatrixIdentity(u2,m_sym,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      cout << "\nTest symmetric multiplications" << endl;
    {
      if (gVerbose)
        cout << "\n  Test m * m_sym == m_sym * m == m_sym * m_sym  multiplications" << endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD mm      = m * m;
      TMatrixD mm_sym1 = m_sym * m_sym;
      TMatrixD mm_sym2 = m     * m_sym;
      TMatrixD mm_sym3 = m_sym * m;
      ok &= VerifyMatrixIdentity(mm,mm_sym1,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mm,mm_sym2,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mm,mm_sym3,gVerbose,epsilon);
    }

    {
      if (gVerbose)
        cout << "\n  Test m^T * m_sym == m_sym^T * m == m_sym^T * m_sym  multiplications" << endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD mtm      = TMatrixD(m    ,TMatrixD::kTransposeMult,m);
      TMatrixD mtm_sym1 = TMatrixD(m_sym,TMatrixD::kTransposeMult,m_sym);
      TMatrixD mtm_sym2 = TMatrixD(m    ,TMatrixD::kTransposeMult,m_sym);
      TMatrixD mtm_sym3 = TMatrixD(m_sym,TMatrixD::kTransposeMult,m);
      ok &= VerifyMatrixIdentity(mtm,mtm_sym1,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mtm,mtm_sym2,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mtm,mtm_sym3,gVerbose,epsilon);
    }

    {
      if (gVerbose)
        cout << "\n  Test m * m_sym^T == m_sym * m^T == m_sym * m_sym^T  multiplications" << endl;
      TMatrixD m = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD mmt      = TMatrixD(m    ,TMatrixD::kMultTranspose,m);
      TMatrixD mmt_sym1 = TMatrixD(m_sym,TMatrixD::kMultTranspose,m_sym);
      TMatrixD mmt_sym2 = TMatrixD(m    ,TMatrixD::kMultTranspose,m_sym);
      TMatrixD mmt_sym3 = TMatrixD(m_sym,TMatrixD::kMultTranspose,m);
      ok &= VerifyMatrixIdentity(mmt,mmt_sym1,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mmt,mmt_sym2,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(mmt,mmt_sym3,gVerbose,epsilon);
    }

    {
      if (gVerbose)
        cout << "\n  Test n * m_sym == n * m multiplications" << endl;
      TMatrixD n = THilbertMatrixD(-1,msize,-1,msize);
      TMatrixD m = n;
      n(1,3) = TMath::Pi();
      n(3,1) = TMath::Pi();
      TMatrixDSym m_sym(-1,msize,m.GetMatrixArray());
      TMatrixD nm1 = n * m_sym;
      TMatrixD nm2 = n * m;
      ok &= VerifyMatrixIdentity(nm1,nm2,gVerbose,epsilon);
    }
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Test inline multiplications by a DiagMatrix" << endl;
    TMatrixDSym ms = THilbertMatrixDSym(msize);
    ms(1,3) = TMath::Pi();
    ms(3,1) = TMath::Pi();
    TVectorD v(msize);
    for (i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = 1+i;
    TMatrixDSym diag(msize);
    TMatrixDDiag d(diag); d = v;
    TMatrixDSym eth = ms;
    for (i = eth.GetRowLwb(); i <= eth.GetRowUpb(); i++)
      for (j = eth.GetColLwb(); j <= eth.GetColUpb(); j++)
        eth(i,j) *= v(j);
    TMatrixD m2 = ms * diag;
    ok &= VerifyMatrixIdentity(m2,eth,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Test XPP = X where P is a permutation matrix" << endl;
    TMatrixDSym ms = THilbertMatrixDSym(msize);
    ms(2,3) = TMath::Pi();
    ms(3,2) = TMath::Pi();
    TMatrixDSym eth = ms;
    TMatrixDSym p(msize);
    for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++)
      p(p.GetRowUpb()+p.GetRowLwb()-i,i) = 1;
    TMatrixD m2 = ms * p;
    m2 *= p;
    ok &= VerifyMatrixIdentity(m2,eth,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Test general matrix multiplication through inline mult" << endl;
    TMatrixDSym ms = THilbertMatrixDSym(msize);
    ms(2,3) = TMath::Pi();
    ms(3,2) = TMath::Pi();
    TMatrixDSym mt(TMatrixDSym::kTransposed,ms);
    TMatrixDSym p = THilbertMatrixDSym(msize);
    TMatrixDDiag(p) += 1;
    TMatrixD mp(ms,TMatrixD::kMult,p);
    TMatrixDSym m1 = ms;
    TMatrixD m3(ms,TMatrixD::kMult,p);
    memcpy(ms.GetMatrixArray(),m3.GetMatrixArray(),msize*msize*sizeof(Double_t));
    ok &= VerifyMatrixIdentity(ms,mp,gVerbose,epsilon);
    TMatrixD mp1(mt,TMatrixD::kTransposeMult,p);
    ok &= VerifyMatrixIdentity(ms,mp1,gVerbose,epsilon);
    ok &= ( !(m1 == ms) ) ? kTRUE : kFALSE;
    TMatrixDSym mp2(TMatrixDSym::kZero,ms);
    ok &= ( mp2 == 0 ) ? kTRUE : kFALSE;

    if (gVerbose)
      cout << "Test XP=X*P  vs XP=X;XP*=P" << endl;
    TMatrixD mp3 = m1*p;
    ok &= VerifyMatrixIdentity(ms,mp3,gVerbose,epsilon);
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Check to see UU' = U'U = E when U is the Haar matrix" << endl;
    {
      const Int_t order = 5;
      const Int_t no_sub_cols = (1<<order)-5;
      TMatrixD haarb = THaarMatrixD(5,no_sub_cols);
      TMatrixD haarb_t(TMatrixD::kTransposed,haarb);
      TMatrixD hth(haarb_t,TMatrixD::kMult,haarb);
      TMatrixDSym  hth1(TMatrixDSym::kAtA,haarb);
      ok &= VerifyMatrixIdentity(hth,hth1,gVerbose,epsilon);
    }

    {
      TMatrixD haar = THaarMatrixD(5);
      TMatrixD unit(TMatrixD::kUnit,haar);
      TMatrixD haar_t(TMatrixD::kTransposed,haar);
      TMatrixDSym hths(TMatrixDSym::kAtA,haar);
      TMatrixD hht(haar,TMatrixD::kMult,haar_t);
      TMatrixD hht1 = haar; hht1 *= haar_t;
      ok &= VerifyMatrixIdentity(unit,hths,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(unit,hht,gVerbose,epsilon);
      ok &= VerifyMatrixIdentity(unit,hht1,gVerbose,epsilon);
    }
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Check to see A.Similarity(Haar) = Haar * A * Haar^T" <<
              " and A.SimilarityT(Haar) = Haar^T * A * Haar" << endl;
    {
      TMatrixD    h  = THaarMatrixD(5);
      TMatrixDSym ms = THilbertMatrixDSym(1<<5);
      TMatrixD    hmht = h*TMatrixD(ms,TMatrixD::kMultTranspose,h);
      ok &= VerifyMatrixIdentity(ms.Similarity(h),hmht,gVerbose,epsilon);
      TMatrixD    htmh = TMatrixD(h,TMatrixD::kTransposeMult,ms)*h;
      ok &= VerifyMatrixIdentity(ms.SimilarityT(h),htmh,gVerbose,epsilon);
    }
    if (gVerbose)
      cout << "Check to see A.Similarity(B_sym) = A.Similarity(B)" << endl;
    {
      TMatrixDSym nsym = THilbertMatrixDSym(5);
      TMatrixD    n    = THilbertMatrixD(5,5);
      TMatrixDSym ms   = THilbertMatrixDSym(5);
      ok &= VerifyMatrixIdentity(ms.Similarity(nsym),ms.Similarity(n),gVerbose,epsilon);
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(11,"Symmetric Matrix Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify vector-matrix multiplications
//
void mstress_vm_multiplications()
{
  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Double_t epsilon = EPSILON*msize;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      cout << "\n---> Verify vector-matrix multiplications "
             "for matrices of the characteristic size " << msize << endl;

    {
      if (verbose)
        cout << "\nCheck shrinking a vector by multiplying by a non-sq unit matrix" << endl;
      TVectorD vb(-2,msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()-i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixD mc(1,msize-2,-2,msize);       // contracting matrix
      mc.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= mc;
      v2.ResizeTo(1,msize-2);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Check expanding a vector by multiplying by a non-sq unit matrix" << endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixD me(2,msize+5,0,msize-1);    // expanding matrix
      me.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= me;
      v2.ResizeTo(v1);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Check general matrix-vector multiplication" << endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      TMatrixD vm(msize,1);
      TMatrixDColumn(vm,0) = vb;
      TMatrixD m = THilbertMatrixD(0,msize,0,msize-1);
      vb *= m;
      ok &= ( vb.GetLwb() == 0 ) ? kTRUE : kFALSE;
      TMatrixD mvm(m,TMatrixD::kMult,vm);
      TMatrixD mvb(msize+1,1);
      TMatrixDColumn(mvb,0) = vb;
      ok &= VerifyMatrixIdentity(mvb,mvm,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Check symmetric matrix-vector multiplication" << endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      TMatrixD vm(msize,1);
      TMatrixDColumn(vm,0) = vb;
      TMatrixDSym ms = THilbertMatrixDSym(0,msize-1);
      vb *= ms;
      ok &= ( vb.GetLwb() == 0 ) ? kTRUE : kFALSE;
      TMatrixD mvm(ms,TMatrixD::kMult,vm);
      TMatrixD mvb(msize,1);
      TMatrixDColumn(mvb,0) = vb;
      ok &= VerifyMatrixIdentity(mvb,mvm,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        cout << "Vector Matrix Multiplications failed for size " << msize << endl;
      break;
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(12,"Matrix Vector Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix inversion
//
void mstress_inversion()
{
  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Double_t epsilon = EPSILON*msize/10;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      cout << "\n---> Verify matrix inversion for square matrices of size " << msize << endl;
    {
      if (verbose)
        cout << "\nTest inversion of a diagonal matrix" << endl;
      TMatrixD m(-1,msize,-1,msize);
      TMatrixD mi(TMatrixD::kZero,m);
      for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
        m(i,i)=i-m.GetRowLwb()+1;
        mi(i,i) = 1/m(i,i);
      }
      TMatrixD mi1(TMatrixD::kInverted,m);
      m.Invert();
      ok &= VerifyMatrixIdentity(m,mi,verbose,epsilon);
      ok &= VerifyMatrixIdentity(mi1,mi,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Test inversion of an orthonormal (Haar) matrix" << endl;
      const Int_t order = Int_t(TMath::Log(msize)/TMath::Log(2));
      TMatrixD mr = THaarMatrixD(order);
      TMatrixD morig = mr;
      TMatrixD mt(TMatrixD::kTransposed,mr);
      Double_t det = -1;         // init to a wrong val to see if it's changed
      mr.Invert(&det);
      ok &= VerifyMatrixIdentity(mr,mt,verbose,epsilon);
      ok &= ( TMath::Abs(det-1) <= msize*epsilon ) ? kTRUE : kFALSE;
      if (verbose) {
        cout << "det = " << det << " deviation= " << TMath::Abs(det-1);
        cout << ( (TMath::Abs(det-1) < msize*epsilon) ? " OK" : " too large") <<endl;
      }
      TMatrixD mti(TMatrixD::kInverted,mt);
      ok &= VerifyMatrixIdentity(mti,morig,verbose,msize*epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Test inversion of a good matrix with diagonal dominance" << endl;
      TMatrixD mh = THilbertMatrixD(msize,msize);
      TMatrixDDiag(mh) += 1;
      TMatrixD morig = mh;
      Double_t det_inv = 0;
      const Double_t det_comp = mh.Determinant();
      mh.Invert(&det_inv);
      if (verbose) {
        cout << "\tcomputed determinant             " << det_comp << endl;
        cout << "\tdeterminant returned by Invert() " << det_inv << endl;
      }

      if (verbose)
        cout << "\tcheck to see M^(-1) * M is E" << endl;
      TMatrixD mim(mh,TMatrixD::kMult,morig);
      TMatrixD unit(TMatrixD::kUnit,mh);
      ok &= VerifyMatrixIdentity(mim,unit,verbose,epsilon);

      if (verbose)
        cout << "Test inversion through the matrix decompositions" << endl;

      TMatrixDSym ms = THilbertMatrixDSym(msize);
      TMatrixDDiag(ms) += 1;
      if (verbose)
        cout << "Test inversion through SVD" << endl;
      TMatrixD inv_svd (msize,msize); TDecompSVD  svd (ms); svd.Invert(inv_svd);
      ok &= VerifyMatrixIdentity(inv_svd,mh,verbose,epsilon);
      if (verbose)
        cout << "Test inversion through LU" << endl;
      TMatrixD inv_lu  (msize,msize); TDecompLU   lu  (ms); lu.Invert(inv_lu);
      ok &= VerifyMatrixIdentity(inv_lu,mh,verbose,epsilon);
      if (verbose)
        cout << "Test inversion through Cholesky" << endl;
      TMatrixDSym inv_chol(msize); TDecompChol chol(ms); chol.Invert(inv_chol);
      ok &= VerifyMatrixIdentity(inv_chol,mh,verbose,epsilon);
      if (verbose)
        cout << "Test inversion through QRH" << endl;
      TMatrixD inv_qrh (msize,msize); TDecompQRH  qrh (ms); qrh.Invert(inv_qrh);
      ok &= VerifyMatrixIdentity(inv_qrh,mh,verbose,epsilon);
      if (verbose)
        cout << "Test inversion through Bunch-Kaufman" << endl;
      TMatrixDSym inv_bk(msize); TDecompBK bk(ms); chol.Invert(inv_bk);
      ok &= VerifyMatrixIdentity(inv_bk,mh,verbose,epsilon);

      if (verbose)
        cout << "\tcheck to see M * M^(-1) is E" << endl;
      TMatrixD mmi = morig; mmi *= mh;
      ok &= VerifyMatrixIdentity(mmi,unit,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        cout << "Matrix Inversion failed for size " << msize << endl;
      break;
    }
  }

  {
    if (gVerbose) {
      cout << "Check to see that Invert() and InvertFast() give identical results" <<endl;
      cout << " for size < (7x7)" <<endl;
    }
    Int_t size;
    for (size = 2; size < 7; size++) {
      TMatrixD m1 = THilbertMatrixD(size,size);
      m1(0,1) = TMath::Pi();
      TMatrixDDiag(m1) += 1;
      TMatrixD m2 = m1;
      Double_t det1 = 0.0;
      Double_t det2 = 0.0;
      m1.Invert(&det1);
      m2.InvertFast(&det2);
      ok &= VerifyMatrixIdentity(m1,m2,gVerbose,EPSILON);
      ok &= (TMath::Abs(det1-det2) < EPSILON);
      if (gVerbose) {
        cout << "det(Invert)= " << det1 << "  det(InvertFast)= " << det2 <<endl;
        cout << " deviation= " << TMath::Abs(det1-det2);
        cout << ( (TMath::Abs(det1-det2) <  EPSILON) ? " OK" : " too large") <<endl;
      }
    }
    for (size = 2; size < 7; size++) {
      TMatrixDSym ms1 = THilbertMatrixDSym(size);
      TMatrixDDiag(ms1) += 1;
      TMatrixDSym ms2 = ms1;
      Double_t det1 = 0.0;
      Double_t det2 = 0.0;
      ms1.Invert(&det1);
      ms2.InvertFast(&det2);
      ok &= VerifyMatrixIdentity(ms1,ms2,gVerbose,EPSILON);
      ok &= (TMath::Abs(det1-det2) < EPSILON);
      if (gVerbose) {
        cout << "det(Invert)= " << det1 << "  det(InvertFast)= " << det2 <<endl;
        cout << " deviation= " << TMath::Abs(det1-det2);
        cout << ( (TMath::Abs(det1-det2) <  EPSILON) ? " OK" : " too large") <<endl;
      }
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(13,"Matrix Inversion",ok);
}

//
//------------------------------------------------------------------------
//           Test matrix I/O
//
void mstress_matrix_io()
{
  if (gVerbose)
    cout << "\n---> Test matrix I/O" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TFile *f = new TFile("stress-vmatrix.root", "RECREATE");

  Char_t name[80];
  Int_t iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m(msize,msize);
    m = pattern;

    Double_t *pattern_array = new Double_t[msize*msize];
    for (Int_t i = 0; i < msize*msize; i++)
      pattern_array[i] = pattern;
    TMatrixD ma;
    ma.Use(msize,msize,pattern_array);

    TMatrixDSym ms(msize);
    ms = pattern;

    if (verbose)
      cout << "\nWrite matrix m to database" << endl;
    snprintf(name,80,"m_%d",msize);
    m.Write(name);

    if (verbose)
      cout << "\nWrite matrix ma which adopts to database" << endl;
    snprintf(name,80,"ma_%d",msize);
    ma.Write(name);

    if (verbose)
      cout << "\nWrite symmetric matrix ms to database" << endl;
    snprintf(name,80,"ms_%d",msize);
    ms.Write(name);

    delete [] pattern_array;

    iloop--;
  }

  if (gVerbose)
    cout << "\nClose database" << endl;
  delete f;

  if (gVerbose)
    cout << "\nOpen database in read-only mode and read matrix" << endl;
  TFile *f1 = new TFile("stress-vmatrix.root");

  iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m(msize,msize);
    m = pattern;
    snprintf(name,80,"m_%d",msize);
    TMatrixD *mr  = (TMatrixD*) f1->Get(name);
    snprintf(name,80,"ma_%d",msize);
    TMatrixD *mar = (TMatrixD*) f1->Get(name);
    snprintf(name,80,"ms_%d",msize);
    TMatrixDSym *msr = (TMatrixDSym*) f1->Get(name);

    if (verbose)
      cout << "\nRead matrix should be same as original" << endl;
    ok &= ((*mr)  == m) ? kTRUE : kFALSE;
    ok &= ((*mar) == m) ? kTRUE : kFALSE;
    ok &= ((*msr) == m) ? kTRUE : kFALSE;

    iloop--;
  }

  delete f1;

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(14,"Matrix Persistence",ok);
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void spstress_allocation(Int_t msize)
{
  if (gVerbose)
    cout << "\n\n---> Test allocation and compatibility check" << endl;

  Int_t i,j;
  Bool_t ok = kTRUE;

  TMatrixDSparse m1(0,3,0,msize-1);
  {
    Int_t nr = 4*msize;
    Int_t    *irow = new Int_t[nr];
    Int_t    *icol = new Int_t[nr];
    Double_t *val  = new Double_t[nr];

    Int_t n = 0;
    for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++) {
      for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++) {
        irow[n] = i;
        icol[n] = j;
        val[n] = TMath::Pi()*i+TMath::E()*j;
        n++;
      }
    }
    m1.SetMatrixArray(nr,irow,icol,val);
    delete [] irow;
    delete [] icol;
    delete [] val;
  }

  TMatrixDSparse m2(0,3,0,msize-1);
  TMatrixDSparse m3(1,4,0,msize-1);
  TMatrixDSparse m4(m1);

  if (gVerbose) {
    cout << "\nStatus information reported for matrix m3:" << endl;
    cout << "  Row lower bound ... " << m3.GetRowLwb() << endl;
    cout << "  Row upper bound ... " << m3.GetRowUpb() << endl;
    cout << "  Col lower bound ... " << m3.GetColLwb() << endl;
    cout << "  Col upper bound ... " << m3.GetColUpb() << endl;
    cout << "  No. rows ..........." << m3.GetNrows()  << endl;
    cout << "  No. cols ..........." << m3.GetNcols()  << endl;
    cout << "  No. of elements ...." << m3.GetNoElements() << endl;
  }

  if (gVerbose)
    cout << "Check matrices 1 & 4 for compatibility" << endl;
  ok &= AreCompatible(m1,m4,gVerbose);

  if (gVerbose)
    cout << "m2 has to be compatible with m3 after resizing to m3" << endl;
  m2.ResizeTo(m3);
  ok &= AreCompatible(m2,m3,gVerbose);

  TMatrixD m5_d(m1.GetNrows()+1,m1.GetNcols()+5);
  for (i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
    for (j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
      m5_d(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixDSparse m5(m5_d);

  if (gVerbose)
    cout << "m1 has to be compatible with m5 after resizing to m5" << endl;
  m1.ResizeTo(m5.GetNrows(),m5.GetNcols());
  ok &= AreCompatible(m1,m5,gVerbose);

  if (gVerbose)
    cout << "m1 has to be equal to m4 after stretching and shrinking" << endl;
  m1.ResizeTo(m4.GetNrows(),m4.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m4,gVerbose,EPSILON);
  if (gVerbose)
    cout << "m5 has to be equal to m1 after shrinking" << endl;
  m5.ResizeTo(m1.GetNrows(),m1.GetNcols());
  ok &= VerifyMatrixIdentity(m1,m5,gVerbose,msize*EPSILON);

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(1,"Allocation, Resizing",ok);
}

//
//------------------------------------------------------------------------
//          Test Filling of matrix
//
void spstress_matrix_fill(Int_t rsize,Int_t csize)
{
  Bool_t ok = kTRUE;

  if (csize < 4) {
    Error("spstress_matrix_fill","rsize should be >= 4");
    ok = kFALSE;
    StatusPrint(2,"Filling, Inserting, Using",ok);
    return;
  }

  if (csize < 4) {
    Error("spstress_matrix_fill","csize should be >= 4");
    ok = kFALSE;
    StatusPrint(2,"Filling, Inserting, Using",ok);
    return;
  }

  TMatrixD m_d(-1,rsize-2,1,csize);
  for (Int_t i = m_d.GetRowLwb(); i <= m_d.GetRowUpb(); i++)
    for (Int_t j = m_d.GetColLwb(); j <= m_d.GetColUpb(); j++)
      // Keep column 2 on zero
      if (j != 2)
        m_d(i,j) = TMath::Pi()*i+TMath::E()*j;
  TMatrixDSparse m(m_d);

  {
    if (gVerbose)
      cout << "Check filling through operator(i,j) without setting sparse index" << endl;
    TMatrixDSparse m1(-1,rsize-2,1,csize);

    for (Int_t i = m1.GetRowLwb(); i <= m1.GetRowUpb(); i++)
      for (Int_t j = m1.GetColLwb(); j <= m1.GetColUpb(); j++)
        if (j != 2)
          m1(i,j) = TMath::Pi()*i+TMath::E()*j;
    ok &= VerifyMatrixIdentity(m1,m,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      cout << "Check filling through operator(i,j)" << endl;
    TMatrixDSparse m2(-1,rsize-2,1,csize);
    m2.SetSparseIndex(m);

    for (Int_t i = m2.GetRowLwb(); i <= m2.GetRowUpb(); i++)
      for (Int_t j = m2.GetColLwb(); j <= m2.GetColUpb(); j++)
        if (j != 2)
          m2(i,j) = TMath::Pi()*i+TMath::E()*j;
    ok &= VerifyMatrixIdentity(m2,m,gVerbose,EPSILON);
  }

  {
    if (gVerbose)
      cout << "Check insertion/extraction of sub-matrices" << endl;
    {
      TMatrixDSparse m_sub1 = m;
      m_sub1.ResizeTo(0,rsize-2,2,csize);
      TMatrixDSparse m_sub2 = m.GetSub(0,rsize-2,2,csize,"");
      ok &= VerifyMatrixIdentity(m_sub1,m_sub2,gVerbose,EPSILON);
    }

    {
      TMatrixDSparse m3(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(-1,rsize-2,1,csize,"");
      m3.SetSub(-1,1,m_part1);
      ok &= VerifyMatrixIdentity(m,m3,gVerbose,EPSILON);
    }

    {
      TMatrixDSparse m4(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(0,rsize-2,2,csize,"");
      TMatrixDSparse m_part2 = m.GetSub(0,rsize-2,1,1,"");
      TMatrixDSparse m_part3 = m.GetSub(-1,-1,2,csize,"");
      TMatrixDSparse m_part4 = m.GetSub(-1,-1,1,1,"");
      m4.SetSub(0,2,m_part1);
      m4.SetSub(0,1,m_part2);
      m4.SetSub(-1,2,m_part3);
      m4.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m4,gVerbose,EPSILON);
    }

    {
      // change the insertion order
      TMatrixDSparse m5(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(0,rsize-2,2,csize,"");
      TMatrixDSparse m_part2 = m.GetSub(0,rsize-2,1,1,"");
      TMatrixDSparse m_part3 = m.GetSub(-1,-1,2,csize,"");
      TMatrixDSparse m_part4 = m.GetSub(-1,-1,1,1,"");
      m5.SetSub(-1,1,m_part4);
      m5.SetSub(-1,2,m_part3);
      m5.SetSub(0,1,m_part2);
      m5.SetSub(0,2,m_part1);
      ok &= VerifyMatrixIdentity(m,m5,gVerbose,EPSILON);
    }

    {
      TMatrixDSparse m6(-1,rsize-2,1,csize);
      TMatrixDSparse m_part1 = m.GetSub(0,rsize-2,2,csize,"S");
      TMatrixDSparse m_part2 = m.GetSub(0,rsize-2,1,1,"S");
      TMatrixDSparse m_part3 = m.GetSub(-1,-1,2,csize,"S");
      TMatrixDSparse m_part4 = m.GetSub(-1,-1,1,1,"S");
      m6.SetSub(0,2,m_part1);
      m6.SetSub(0,1,m_part2);
      m6.SetSub(-1,2,m_part3);
      m6.SetSub(-1,1,m_part4);
      ok &= VerifyMatrixIdentity(m,m6,gVerbose,EPSILON);
    }
  }

  {
    if (gVerbose)
      cout << "Check insertion/extraction of rows" << endl;

    TMatrixDSparse m1 = m;
    TVectorD v1(1,csize);
    TVectorD v2(1,csize);
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      v1 = TMatrixDSparseRow(m,i);
      m1.InsertRow(i,1,v1.GetMatrixArray());
      ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
      m1.InsertRow(i,3,v1.GetMatrixArray()+2,v1.GetNrows()-2);
      ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

      m1.ExtractRow(i,1,v2.GetMatrixArray());
      ok &= VerifyVectorIdentity(v1,v2,gVerbose,EPSILON);
      m1.ExtractRow(i,3,v2.GetMatrixArray()+2,v1.GetNrows()-2);
      ok &= VerifyVectorIdentity(v1,v2,gVerbose,EPSILON);
    }
  }

  {
    if (gVerbose)
      cout << "Check array Use" << endl;
    {
      TMatrixDSparse *m1a = new TMatrixDSparse(m);
      TMatrixDSparse *m2a = new TMatrixDSparse();
      m2a->Use(*m1a);
      m2a->Sqr();
      TMatrixDSparse m7 = m; m7.Sqr();
      ok &= VerifyMatrixIdentity(m7,*m1a,gVerbose,EPSILON);
      delete m1a;
      delete m2a;
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(2,"Filling, Inserting, Using",ok);
}

//
//------------------------------------------------------------------------
//                Test uniform element operations
//
void spstress_element_op(Int_t rsize,Int_t csize)
{
  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TMatrixDSparse m(-1,rsize-2,1,csize);

  if (gVerbose)
    cout << "Creating zero m1 ..." << endl;
  TMatrixDSparse m1(TMatrixDSparse::kZero, m);
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Comparing m1 with 0 ..." << endl;
  R__ASSERT(m1 == 0);
  R__ASSERT(!(m1 != 0));

  if (gVerbose)
    cout << "Writing a pattern " << pattern << " by assigning through SetMatrixArray..." << endl;
  {
    const Int_t nr = rsize*csize;
    Int_t    *irow = new Int_t[nr];
    Int_t    *icol = new Int_t[nr];
    Double_t *val  = new Double_t[nr];

    Int_t n = 0;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        irow[n] = i;
        icol[n] = j;
        val[n] = pattern;
        n++;
      }
    }
    m.SetMatrixArray(nr,irow,icol,val);
    delete [] irow;
    delete [] icol;
    delete [] val;
    ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "Writing the pattern by assigning to m1 as a whole ..."  << endl;
  m1.SetSparseIndex(m);
  m1 = pattern;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Comparing m and m1 ..." << endl;
  R__ASSERT(m == m1);
  if (gVerbose)
    cout << "Comparing (m=0) and m1 ..." << endl;
  R__ASSERT(!((m=0) == m1));

  if (gVerbose)
    cout << "Clearing m1 ..." << endl;
  m1.Zero();
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nSet m = pattern" << endl;
  m = pattern;
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   add the doubled pattern with the negative sign" << endl;
  m += -2*pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   subtract the trippled pattern with the negative sign" << endl;
  m -= -3*pattern;
  ok &= VerifyMatrixValue(m,2*pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nVerify comparison operations when all elems are the same" << endl;
  m = pattern;
  R__ASSERT( m == pattern && !(m != pattern) );
  R__ASSERT( m > 0 && m >= pattern && m <= pattern );
  R__ASSERT( m > -pattern && m >= -pattern );
  R__ASSERT( m <= pattern && !(m < pattern) );
  m -= 2*pattern;
  R__ASSERT( m  < -pattern/2 && m <= -pattern/2 );
  R__ASSERT( m  >= -pattern && !(m > -pattern) );

  if (gVerbose)
    cout << "\nVerify comparison operations when not all elems are the same" << endl;
  {
    Int_t nr = rsize*csize;
    Int_t    *irow = new Int_t[nr];
    Int_t    *icol = new Int_t[nr];
    Double_t *val  = new Double_t[nr];

    Int_t n = 0;
    for (Int_t i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
      for (Int_t j = m.GetColLwb(); j <= m.GetColUpb(); j++) {
        irow[n] = i;
        icol[n] = j;
        val[n] = pattern;
        n++;
      }
    }
    val[n-1] = pattern-1;
    m.SetMatrixArray(nr,irow,icol,val);
    delete [] irow;
    delete [] icol;
    delete [] val;
  }

  R__ASSERT( !(m == pattern) && !(m != pattern) );
  R__ASSERT( m != 0 );                   // none of elements are 0
  R__ASSERT( !(m >= pattern) && m <= pattern && !(m<pattern) );
  R__ASSERT( !(m <= pattern-1) && m >= pattern-1 && !(m>pattern-1) );
  if (gVerbose)
    cout << "\nAssign 2*pattern to m by repeating additions" << endl;
  m = 0; m += pattern; m += pattern;
  if (gVerbose)
    cout << "Assign 2*pattern to m1 by multiplying by two " << endl;
  m1.SetSparseIndex(m);
  m1 = pattern; m1 *= 2;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    cout << "Multiply m1 by one half returning it to the 1*pattern" << endl;
  m1 *= 1/2.;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nAssign -pattern to m and m1" << endl;
  m = 0; m -= pattern; m1 = -pattern;
  ok &= VerifyMatrixValue(m,-pattern,gVerbose,EPSILON);
  R__ASSERT( m == m1 );
  if (gVerbose)
    cout << "m = sqrt(sqr(m)); m1 = abs(m1); Now m and m1 have to be the same" << endl;
  m.Sqr();
  ok &= VerifyMatrixValue(m,pattern*pattern,gVerbose,EPSILON);
  m.Sqrt();
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  m1.Abs();
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  ok &= VerifyMatrixIdentity(m1,m,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(3,"Uniform matrix operations",ok);
}

//
//------------------------------------------------------------------------
//        Test binary matrix element-by-element operations
//
void spstress_binary_ebe_op(Int_t rsize, Int_t csize)
{
  if (gVerbose)
    cout << "\n---> Test Binary Matrix element-by-element operations" << endl;

  Bool_t ok = kTRUE;

  const Double_t pattern = 4.25;

  TMatrixD m_d(2,rsize+1,0,csize-1); m_d = 1;
  TMatrixDSparse m (TMatrixDSparse::kZero,m_d); m.SetSparseIndex (m_d);
  TMatrixDSparse m1(TMatrixDSparse::kZero,m);   m1.SetSparseIndex(m_d);

  TMatrixDSparse mp(TMatrixDSparse::kZero,m1);  mp.SetSparseIndex(m_d);
  {
    for (Int_t i = mp.GetRowLwb(); i <= mp.GetRowUpb(); i++)
      for (Int_t j = mp.GetColLwb(); j <= mp.GetColUpb(); j++)
        mp(i,j) = TMath::Pi()*i+TMath::E()*(j+1);
  }

  if (gVerbose)
    cout << "\nVerify assignment of a matrix to the matrix" << endl;
  m  = pattern;
  m1 = m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  R__ASSERT( m1 == m );

  if (gVerbose)
    cout << "\nAdding the matrix to itself, uniform pattern " << pattern << endl;
  m.Zero(); m.SetSparseIndex(m_d); m = pattern;

  m1 = m; m1 += m1;
  ok &= VerifyMatrixValue(m1,2*pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  subtracting two matrices ..." << endl;
  m1 -= m;
  ok &= VerifyMatrixValue(m1,pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  subtracting the matrix from itself" << endl;
  m1 -= m1;
  ok &= VerifyMatrixValue(m1,0.,gVerbose,EPSILON);
  m1.SetSparseIndex(m_d);

  if (gVerbose) {
    cout << "\nArithmetic operations on matrices with not the same elements" << endl;
    cout << "   adding mp to the zero matrix..." << endl;
  }
  m.Zero(); m.SetSparseIndex(m_d); m += mp;
  ok &= VerifyMatrixIdentity(m,mp,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   making m = 3*mp and m1 = 3*mp, via add() and succesive mult" << endl;
  m1 = m;
  Add(m,2.,mp);
  m1 += m1; m1 += mp;
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   clear both m and m1, by subtracting from itself and via add()" << endl;
  m1 -= m1;
  Add(m,-3.,mp);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);

  if (gVerbose) {
    cout << "\nTesting element-by-element multiplications and divisions" << endl;
    cout << "   squaring each element with sqr() and via multiplication" << endl;
  }
  m.SetSparseIndex(m_d);  m = mp;
  m1.SetSparseIndex(m_d); m1 = mp;
  m.Sqr();
  ElementMult(m1,m1);
  ok &= VerifyMatrixIdentity(m,m1,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   compare (m = pattern^2)/pattern with pattern" << endl;
  m = pattern; m1 = pattern;
  m.Sqr();
  ElementDiv(m,m1);
  ok &= VerifyMatrixValue(m,pattern,gVerbose,EPSILON);
  if (gVerbose)
    Compare(m1,m);

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(4,"Binary Matrix element-by-element operations",ok);
}

//
//------------------------------------------------------------------------
//              Verify matrix transposition
//
void spstress_transposition(Int_t msize)
{
  if (gVerbose) {
    cout << "\n---> Verify matrix transpose "
            "for matrices of a characteristic size " << msize << endl;
  }

  Bool_t ok = kTRUE;
  {
    if (gVerbose)
      cout << "\nCheck to see that a square UnitMatrix stays the same";
    TMatrixDSparse m(msize,msize);
    m.UnitMatrix();
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "\nTest a non-square UnitMatrix";
    TMatrixDSparse m(msize,msize+1);
    m.UnitMatrix();
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows() );
    const Int_t rowlwb = m.GetRowLwb();
    const Int_t collwb = m.GetColLwb();
    const Int_t upb = TMath::Min(m.GetRowUpb(),m.GetColUpb());
    TMatrixDSparse m_part  = m.GetSub(rowlwb,upb,collwb,upb);
    TMatrixDSparse mt_part = mt.GetSub(rowlwb,upb,collwb,upb);
    ok &= ( m_part == mt_part ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "\nCheck to see that a symmetric (Hilbert)Matrix stays the same";
    TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize,msize));
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    ok &= ( m == mt ) ? kTRUE : kFALSE;
  }

  {
    if (gVerbose)
      cout << "\nCheck transposing a non-symmetric matrix";
    TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize+1,msize));
    m(1,2) = TMath::Pi();
    TMatrixDSparse mt(TMatrixDSparse::kTransposed,m);
    R__ASSERT(m.GetNrows() == mt.GetNcols() && m.GetNcols() == mt.GetNrows());
    R__ASSERT(mt(2,1)  == (Double_t)TMath::Pi() && mt(1,2)  != (Double_t)TMath::Pi());
    R__ASSERT(mt[2][1] == (Double_t)TMath::Pi() && mt[1][2] != (Double_t)TMath::Pi());

    if (gVerbose)
      cout << "\nCheck double transposing a non-symmetric matrix" << endl;
    TMatrixDSparse mtt(TMatrixDSparse::kTransposed,mt);
    ok &= ( m == mtt ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(5,"Matrix transposition",ok);
}

//
//------------------------------------------------------------------------
//             Verify the norm calculation
//
void spstress_norms(Int_t rsize_req,Int_t csize_req)
{
  if (gVerbose)
    cout << "\n---> Verify norm calculations" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 10.25;

  Int_t rsize = rsize_req;
  if (rsize%2 != 0)  rsize--;
  Int_t csize = csize_req;
  if (csize%2 != 0)  csize--;
  if (rsize%2 == 1 || csize%2 == 1) {
    cout << "rsize: " << rsize <<endl;
    cout << "csize: " << csize <<endl;
    Fatal("spstress_norms","Sorry, size of the matrix to test must be even for this test\n");
  }

  TMatrixD m_d(rsize,csize); m_d = 1;
  TMatrixDSparse m(rsize,csize); m.SetSparseIndex(m_d);

  if (gVerbose)
    cout << "\nAssign " << pattern << " to all the elements and check norms" << endl;
  m = pattern;
  if (gVerbose)
    cout << "  1. (col) norm should be pattern*nrows" << endl;
  ok &= ( m.Norm1() == pattern*m.GetNrows() ) ? kTRUE : kFALSE;
  ok &= ( m.Norm1() == m.ColNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Inf (row) norm should be pattern*ncols" << endl;
  ok &= ( m.NormInf() == pattern*m.GetNcols() ) ? kTRUE : kFALSE;
  ok &= ( m.NormInf() == m.RowNorm() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Square of the Eucl norm has got to be pattern^2 * no_elems" << endl;
  ok &= ( m.E2Norm() == (pattern*pattern)*m.GetNoElements() ) ? kTRUE : kFALSE;
  TMatrixDSparse m1(m); m1 = 1;
  ok &= ( m.E2Norm() == E2Norm(m+1.,m1) ) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(6,"Matrix Norms",ok);
}

//
//------------------------------------------------------------------------
//               Verify matrix multiplications
//
void spstress_mm_multiplications()
{
  Bool_t ok = kTRUE;
  Int_t i;

  Int_t iloop = 0;
  Int_t nr    = 0;
  while (iloop <= gNrLoop) {
    const Int_t msize = gSizeA[iloop];
    const Double_t epsilon = EPSILON*msize/100;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (msize <= 5) {
       iloop++;
       continue;  // some references to m(3,..)
    }

    if (verbose)
      cout << "\n---> Verify matrix multiplications "
              "for matrices of the characteristic size " << msize << endl;

    {
      if (verbose)
        cout << "\nTest inline multiplications of the UnitMatrix" << endl;
      if (ok)
      {
        TMatrixDSparse m = TMatrixD(THilbertMatrixD(-1,msize,-1,msize));
        TMatrixDSparse ur(TMatrixDSparse::kUnit,m);
        TMatrixDSparse ul(TMatrixDSparse::kUnit,m);
        m(3,1) = TMath::Pi();
        ul *= m;
        m  *= ur;
        ok &= VerifyMatrixIdentity(ul,m,verbose,epsilon);
      }

      if (ok)
      {
        TMatrixD m_d = THilbertMatrixD(-1,msize,-1,msize);
        TMatrixDSparse ur(TMatrixDSparse::kUnit,m_d);
        TMatrixDSparse ul(TMatrixDSparse::kUnit,m_d);
        m_d(3,1) = TMath::Pi();
        ul *= m_d;
        m_d  *= TMatrixD(ur);
        ok &= VerifyMatrixIdentity(ul,m_d,verbose,epsilon);
      }
    }

    if (ok)
    {
      if (verbose)
        cout << "Test XPP = X where P is a permutation matrix" << endl;
      TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize-1,msize));
      m(2,3) = TMath::Pi();
      TMatrixDSparse eth = m;
      TMatrixDSparse p(msize,msize);
      {
        Int_t    *irow = new Int_t[msize];
        Int_t    *icol = new Int_t[msize];
        Double_t *val  = new Double_t[msize];

        Int_t n = 0;
        for (i = p.GetRowLwb(); i <= p.GetRowUpb(); i++) {
          irow[n] = p.GetRowUpb()+p.GetRowLwb()-i;
          icol[n] = i;
          val[n] = 1;
          n++;
        }
        p.SetMatrixArray(msize,irow,icol,val);
        delete [] irow;
        delete [] icol;
        delete [] val;
      }
      m *= p;
      m *= p;
      ok &= VerifyMatrixIdentity(m,eth,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Test general matrix multiplication through inline mult" << endl;
      TMatrixDSparse m = TMatrixD(THilbertMatrixD(msize-2,msize));
      m(3,3) = TMath::Pi();
      TMatrixD p_d = THilbertMatrixD(msize,msize);
      TMatrixDDiag(p_d) += 1;
      TMatrixDSparse mp(m,TMatrixDSparse::kMult,p_d);
      TMatrixDSparse m1 = m;
      m *= p_d;
      ok &= VerifyMatrixIdentity(m,mp,verbose,epsilon);
      TMatrixDSparse pt_d(TMatrixDSparse::kTransposed,p_d);
      TMatrixDSparse mp1(m1,TMatrixDSparse::kMultTranspose,pt_d);
      VerifyMatrixIdentity(m,mp1,verbose,epsilon);

      ok &= ( !(m1 == m) );
      TMatrixDSparse mp2(TMatrixDSparse::kZero,m1);
      ok &= ( mp2 == 0 );
      mp2.SetSparseIndex(m1);
      mp2.Mult(m1,p_d);
      ok &= VerifyMatrixIdentity(m,mp2,verbose,epsilon);

      if (verbose)
        cout << "Test XP=X*P  vs XP=X;XP*=P" << endl;
      TMatrixDSparse mp3 = m1*p_d;
      ok &= VerifyMatrixIdentity(m,mp3,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop++;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        cout << "General Matrix Multiplications failed for size " << msize << endl;
      break;
    }
  }

  if (ok)
  {
    if (gVerbose)
      cout << "Check to see UU' = U'U = E when U is the Haar matrix" << endl;
    const Int_t order = 5;
    const Int_t no_sub_cols = (1<<order)-order;
    TMatrixDSparse haar_sub = TMatrixD(THaarMatrixD(order,no_sub_cols));
    TMatrixDSparse haar_sub_t(TMatrixDSparse::kTransposed,haar_sub);
    TMatrixDSparse hsths(haar_sub_t,TMatrixDSparse::kMult,haar_sub);

    TMatrixDSparse square(no_sub_cols,no_sub_cols);
    TMatrixDSparse units(TMatrixDSparse::kUnit,square);
    ok &= ( hsths.GetNrows() == no_sub_cols && hsths.GetNcols() == no_sub_cols );
    ok &= VerifyMatrixIdentity(hsths,units,gVerbose,EPSILON);

    TMatrixDSparse haar = TMatrixD(THaarMatrixD(order));
    TMatrixDSparse haar_t(TMatrixDSparse::kTransposed,haar);
    TMatrixDSparse unit(TMatrixDSparse::kUnit,haar);

    TMatrixDSparse hth(haar_t,TMatrixDSparse::kMult,haar);
    ok &= VerifyMatrixIdentity(hth,unit,gVerbose,EPSILON);

    TMatrixDSparse hht(haar,TMatrixDSparse::kMultTranspose,haar);
    ok &= VerifyMatrixIdentity(hht,unit,gVerbose,EPSILON);

    TMatrixDSparse hht1 = haar; hht1 *= haar_t;
    ok &= VerifyMatrixIdentity(hht1,unit,gVerbose,EPSILON);

    TMatrixDSparse hht2(TMatrixDSparse::kZero,haar);
    hht2.SetSparseIndex(hht1);
    hht2.Mult(haar,haar_t);
    ok &= VerifyMatrixIdentity(hht2,unit,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(7,"General Matrix Multiplications",ok);
}

//
//------------------------------------------------------------------------
//               Verify vector-matrix multiplications
//
void spstress_vm_multiplications()
{
  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Double_t epsilon = EPSILON*msize;

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      cout << "\n---> Verify vector-matrix multiplications "
             "for matrices of the characteristic size " << msize << endl;

    {
      if (verbose)
        cout << "\nCheck shrinking a vector by multiplying by a non-sq unit matrix" << endl;
      TVectorD vb(-2,msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()-i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixDSparse mc(1,msize-2,-2,msize);       // contracting matrix
      mc.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= mc;
      v2.ResizeTo(1,msize-2);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Check expanding a vector by multiplying by a non-sq unit matrix" << endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      ok &= ( vb != 0 ) ? kTRUE : kFALSE;
      TMatrixDSparse me(2,msize+5,0,msize-1);    // expanding matrix
      me.UnitMatrix();
      TVectorD v1 = vb;
      TVectorD v2 = vb;
      v1 *= me;
      v2.ResizeTo(v1);
      ok &= VerifyVectorIdentity(v1,v2,verbose,epsilon);
    }

    if (ok)
    {
      if (verbose)
        cout << "Check general matrix-vector multiplication" << endl;
      TVectorD vb(msize);
      for (Int_t i = vb.GetLwb(); i <= vb.GetUpb(); i++)
        vb(i) = TMath::Pi()+i;
      TMatrixD vm(msize,1);
      TMatrixDColumn(vm,0) = vb;
      TMatrixD hilbert_with_zeros = THilbertMatrixD(0,msize,0,msize-1);
      TMatrixDRow   (hilbert_with_zeros,3) = 0.0;
      TMatrixDColumn(hilbert_with_zeros,3) = 0.0;
      const TMatrixDSparse m = hilbert_with_zeros;
      vb *= m;
      ok &= ( vb.GetLwb() == 0 ) ? kTRUE : kFALSE;
      TMatrixDSparse mvm(m,TMatrixDSparse::kMult,vm);
      TMatrixD mvb(msize+1,1);
      TMatrixDColumn(mvb,0) = vb;
      ok &= VerifyMatrixIdentity(mvb,mvm,verbose,epsilon);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        cout << "Vector Matrix Multiplications failed for size " << msize << endl;
      break;
    }
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(8,"Matrix Vector Multiplications",ok);
}

//
//------------------------------------------------------------------------
//           Test operations with vectors and sparse matrix slices
//
void spstress_matrix_slices(Int_t vsize)
{
  if (gVerbose)
    cout << "\n---> Test operations with vectors and sparse matrix slices" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TVectorD vc(0,vsize);
  TVectorD vr(0,vsize+1);
  TMatrixD       m_d(0,vsize,0,vsize+1); m_d = pattern;
  TMatrixDSparse m(m_d);

  Int_t i,j;
  if (gVerbose)
    cout << "\nCheck modifying the matrix row-by-row" << endl;
  m = pattern;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    TMatrixDSparseRow(m,i) = pattern+2;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    vr = TMatrixDSparseRow(m,i);
    ok &= VerifyVectorValue(vr,pattern+2,gVerbose,EPSILON);
    vr = TMatrixDSparseRow(m,i+1 > m.GetRowUpb() ? m.GetRowLwb() : i+1);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    TMatrixDSparseRow(m,i) += -2;
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    vr = pattern-2;
    TMatrixDSparseRow(m,i) = vr;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    {
      TMatrixDSparseRow mr(m,i);
      for (j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        mr(j) *= 8;
    }
    vr = TMatrixDSparseRow(m,i);
    ok &= VerifyVectorValue(vr,8*(pattern-2),gVerbose,EPSILON);
    TMatrixDSparseRow(m,i) *= 1./8;
    TMatrixDSparseRow(m,i) += 2;
    vr = TMatrixDSparseRow(m,i);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nCheck modifying the matrix diagonal" << endl;
  m = pattern;
  TMatrixDSparseDiag td = m;
  td = pattern-3;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  vc = TMatrixDSparseDiag(m);
  ok &= VerifyVectorValue(vc,pattern-3,gVerbose,EPSILON);
  td += 3;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  vc = pattern+3;
  td = vc;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  {
    TMatrixDSparseDiag md(m);
    for (j = 0; j < md.GetNdiags(); j++)
      md(j) /= 1.5;
  }
  vc = TMatrixDSparseDiag(m);
  ok &= VerifyVectorValue(vc,(pattern+3)/1.5,gVerbose,EPSILON);
  TMatrixDSparseDiag(m) *= 1.5;
  TMatrixDSparseDiag(m) += -3;
  vc = TMatrixDSparseDiag(m);
  ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
  ok &= ( m == pattern ) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(9,"Matrix Slices to Vectors",ok);
}

//
//------------------------------------------------------------------------
//           Test matrix I/O
//
void spstress_matrix_io()
{
  if (gVerbose)
    cout << "\n---> Test matrix I/O" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TFile *f = new TFile("stress-vmatrix.root", "RECREATE");

  Char_t name[80];
  Int_t iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m_d(msize,msize); m_d = pattern;
    TMatrixDSparse m(m_d);
    TMatrixDSparse ma;
    ma.Use(m);

    if (verbose)
      cout << "\nWrite matrix m to database" << endl;
    snprintf(name,80,"m_%d",msize);
    m.Write(name);

    if (verbose)
      cout << "\nWrite matrix ma which adopts to database" << endl;
    snprintf(name,80,"ma_%d",msize);
    ma.Write(name);

    iloop--;
  }

  if (gVerbose)
    cout << "\nClose database" << endl;
  delete f;

  if (gVerbose)
    cout << "\nOpen database in read-only mode and read matrix" << endl;
  TFile *f1 = new TFile("stress-vmatrix.root");

  iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TMatrixD m_d(msize,msize); m_d = pattern;
    TMatrixDSparse m(m_d);
    snprintf(name,80,"m_%d",msize);
    TMatrixDSparse *mr  = (TMatrixDSparse*) f1->Get(name);
    snprintf(name,80,"ma_%d",msize);
    TMatrixDSparse *mar = (TMatrixDSparse*) f1->Get(name);
    snprintf(name,80,"ms_%d",msize);

    if (verbose)
      cout << "\nRead matrix should be same as original" << endl;
    ok &= ((*mr)  == m) ? kTRUE : kFALSE;
    ok &= ((*mar) == m) ? kTRUE : kFALSE;

    iloop--;
  }

  delete f1;

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(10,"Matrix Persistence",ok);
}

//------------------------------------------------------------------------
//          Test allocation functions and compatibility check
//
void vstress_allocation(Int_t msize)
{
  if (gVerbose)
    cout << "\n\n---> Test allocation and compatibility check" << endl;

  Int_t i;
  Bool_t ok = kTRUE;
  TVectorD v1(msize);
  TVectorD v2(0,msize-1);
  TVectorD v3(1,msize);
  TVectorD v4(v1);

  if (gVerbose) {
    cout << "\nStatus information reported for vector v3:" << endl;
    cout << "  Lower bound ... " << v3.GetLwb() << endl;
    cout << "  Upper bound ... " << v3.GetUpb() << endl;
    cout << "  No. of elements " << v3.GetNoElements() << endl;
  }

  if (gVerbose)
    cout << "\nCheck vectors 1 & 2 for compatibility" << endl;
  ok &= AreCompatible(v1,v2,gVerbose);

  if (gVerbose)
    cout << "Check vectors 1 & 4 for compatibility" << endl;
  ok &= AreCompatible(v1,v4,gVerbose);

  if (gVerbose)
    cout << "v2 has to be compatible with v3 after resizing to v3" << endl;
  v2.ResizeTo(v3);
  ok &= AreCompatible(v2,v3,gVerbose);

  TVectorD v5(v1.GetUpb()+5);
  if (gVerbose)
    cout << "v1 has to be compatible with v5 after resizing to v5.upb" << endl;
  v1.ResizeTo(v5.GetNoElements());
  ok &= AreCompatible(v1,v5,gVerbose);

  {
    if (gVerbose)
      cout << "Check that shrinking does not change remaining elements" << endl;
    TVectorD vb(-1,msize);
    for (i = vb.GetLwb(); i <= vb.GetUpb(); i++)
      vb(i) = i+TMath::Pi();
    TVectorD v = vb;
    ok &= ( v == vb ) ? kTRUE : kFALSE;
    ok &= ( v != 0 ) ? kTRUE : kFALSE;
    v.ResizeTo(0,msize/2);
    for (i = v.GetLwb(); i <= v.GetUpb(); i++)
      ok &= ( v(i) == vb(i) ) ? kTRUE : kFALSE;
    if (gVerbose)
      cout << "Check that expansion expands by zeros" << endl;
    const Int_t old_nelems = v.GetNoElements();
    const Int_t old_lwb    = v.GetLwb();
    v.ResizeTo(vb);
    ok &= ( !(v == vb) ) ? kTRUE : kFALSE;
    for (i = old_lwb; i < old_lwb+old_nelems; i++)
      ok &= ( v(i) == vb(i) ) ? kTRUE : kFALSE;
    for (i = v.GetLwb(); i < old_lwb; i++)
      ok &= ( v(i) == 0 ) ? kTRUE : kFALSE;
    for (i = old_lwb+old_nelems; i <= v.GetUpb(); i++)
      ok &= ( v(i) == 0 ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(1,"Allocation, Filling, Resizing",ok);
}

//
//------------------------------------------------------------------------
//                Test uniform element operations
//
class SinAction : public TElementActionD {
  void Operation(Double_t &element) const { element = TMath::Sin(element); }
  public:
    SinAction() { }
};

class CosAction : public TElementPosActionD {
  Double_t factor;
  void Operation(Double_t &element) const { element = TMath::Cos(factor*fI); }
  public:
    CosAction(Int_t no_elems): factor(2*TMath::Pi()/no_elems) { }
};

void vstress_element_op(Int_t vsize)
{
  if (gVerbose)
    cout << "\n---> Test operations that treat each element uniformly" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TVectorD v(-1,vsize-2);
  TVectorD v1(v);

  if (gVerbose)
    cout << "\nWriting zeros to v..." << endl;
  for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
    v(i) = 0;
  ok &= VerifyVectorValue(v,0.0,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Clearing v1 ..." << endl;
  v1.Zero();
  ok &= VerifyVectorValue(v1,0.0,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Comparing v1 with 0 ..." << endl;
  ok &= (v1 == 0) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "Writing a pattern " << pattern << " by assigning to v(i)..." << endl;
  {
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = pattern;
    ok &= VerifyVectorValue(v,pattern,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "Writing the pattern by assigning to v1 as a whole ..." << endl;
  v1 = pattern;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "Comparing v and v1 ..." << endl;
  ok &= (v == v1) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "Comparing (v=0) and v1 ..." << endl;
  ok &= (!(v.Zero() == v1)) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nClear v and add the pattern" << endl;
  v.Zero();
  v += pattern;
  ok &= VerifyVectorValue(v,pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   add the doubled pattern with the negative sign" << endl;
  v += -2*pattern;
  ok &= VerifyVectorValue(v,-pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "   subtract the trippled pattern with the negative sign" << endl;
  v -= -3*pattern;
  ok &= VerifyVectorValue(v,2*pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nVerify comparison operations" << endl;
  v = pattern;
  ok &= ( v == pattern && !(v != pattern) && v >= pattern && v <= pattern ) ? kTRUE : kFALSE;
  ok &= ( v > 0 && v >= 0 ) ? kTRUE : kFALSE;
  ok &= ( v > -pattern && v >= -pattern ) ? kTRUE : kFALSE;
  ok &= ( v < pattern+1 && v <= pattern+1 ) ? kTRUE : kFALSE;
  v(v.GetUpb()) += 1;
  ok &= ( !(v==pattern)      && !(v != pattern)  && v != pattern-1 ) ? kTRUE : kFALSE;
  ok &= ( v >= pattern       && !(v > pattern)   && !(v >= pattern+1) ) ? kTRUE : kFALSE;
  ok &= ( v <= pattern+1.001 && !(v < pattern+1) && !(v <= pattern) ) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nAssign 2*pattern to v by repeating additions" << endl;
  v = 0; v += pattern; v += pattern;
  if (gVerbose)
    cout << "Assign 2*pattern to v1 by multiplying by two" << endl;
  v1 = pattern; v1 *= 2;
  ok &= VerifyVectorValue(v1,2*pattern,gVerbose,EPSILON);
  ok &= ( v == v1 ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "Multiply v1 by one half returning it to the 1*pattern" << endl;
  v1 *= 1/2.;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nAssign -pattern to v and v1" << endl;
  v.Zero(); v -= pattern; v1 = -pattern;
  ok &= VerifyVectorValue(v,-pattern,gVerbose,EPSILON);
  ok &= ( v == v1 ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "v = sqrt(sqr(v)); v1 = abs(v1); Now v and v1 have to be the same" << endl;
  v.Sqr();
  ok &= VerifyVectorValue(v,pattern*pattern,gVerbose,EPSILON);
  v.Sqrt();
  ok &= VerifyVectorValue(v,pattern,gVerbose,EPSILON);
  v1.Abs();
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);
  ok &= ( v == v1 ) ? kTRUE : kFALSE;

  {
    if (gVerbose)
      cout << "\nCheck out to see that sin^2(x) + cos^2(x) = 1" << endl;
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = 2*TMath::Pi()/v.GetNoElements()*i;
#ifndef __CINT__
    SinAction s;
    v.Apply(s);
    CosAction c(v.GetNoElements());
    v1.Apply(c);
#else
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++) {
      v(i)  = TMath::Sin(v(i));
      v1(i) = TMath::Cos(2*TMath::Pi()/v1.GetNrows()*i);
    }
#endif
    TVectorD v2 = v;
    TVectorD v3 = v1;
    v.Sqr();
    v1.Sqr();
    v += v1;
    ok &= VerifyVectorValue(v,1.,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "\nVerify constructor with initialization" << endl;
#ifndef __CINT__
  TVectorD vi(0,4,0.0,1.0,2.0,3.0,4.0,"END");
#else
  Double_t vval[] = {0.0,1.0,2.0,3.0,4.0};
  TVectorD vi(5,vval);
#endif
  TVectorD vit(5);
  {
    for (Int_t i = vit.GetLwb(); i <= vit.GetUpb(); i++)
      vit(i) = Double_t(i);
    ok &= VerifyVectorIdentity(vi,vit,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(2,"Uniform vector operations",ok);
}

//
//------------------------------------------------------------------------
//                 Test binary vector operations
//
void vstress_binary_op(Int_t vsize)
{
  if (gVerbose)
    cout << "\n---> Test Binary Vector operations" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  const Double_t epsilon = EPSILON*vsize/10;

  TVectorD v(2,vsize+1);
  TVectorD v1(v);

  if (gVerbose)
    cout << "\nVerify assignment of a vector to the vector" << endl;
  v = pattern;
  v1.Zero();
  v1 = v;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);
  ok &= ( v1 == v ) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nAdding one vector to itself, uniform pattern " << pattern << endl;
  v.Zero(); v = pattern;
  v1 = v; v1 += v1;
  ok &= VerifyVectorValue(v1,2*pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  subtracting two vectors ..." << endl;
  v1 -= v;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  subtracting the vector from itself" << endl;
  v1 -= v1;
  ok &= VerifyVectorValue(v1,0.,gVerbose,EPSILON);
  if (gVerbose)
    cout << "  adding two vectors together" << endl;
  v1 += v;
  ok &= VerifyVectorValue(v1,pattern,gVerbose,EPSILON);

  TVectorD vp(2,vsize+1);
  {
    for (Int_t i = vp.GetLwb(); i <= vp.GetUpb(); i++)
      vp(i) = (i-vp.GetNoElements()/2.)*pattern;
  }

  if (gVerbose) {
    cout << "\nArithmetic operations on vectors with not the same elements" << endl;
    cout << "   adding vp to the zero vector..." << endl;
  }
  v.Zero();
  ok &= ( v == 0.0 ) ? kTRUE : kFALSE;
  v += vp;
  ok &= VerifyVectorIdentity(v,vp,gVerbose,epsilon);
  v1 = v;
  if (gVerbose)
    cout << "   making v = 3*vp and v1 = 3*vp, via add() and succesive mult" << endl;
  Add(v,2.,vp);
  v1 += v1; v1 += vp;
  ok &= VerifyVectorIdentity(v,v1,gVerbose,epsilon);
  if (gVerbose)
    cout << "   clear both v and v1, by subtracting from itself and via add()" << endl;
  v1 -= v1;
  Add(v,-3.,vp);
  ok &= VerifyVectorIdentity(v,v1,gVerbose,epsilon);

  if (gVerbose) {
    cout << "\nTesting element-by-element multiplications and divisions" << endl;
    cout << "   squaring each element with sqr() and via multiplication" << endl;
  }
  v = vp; v1 = vp;
  v.Sqr();
  ElementMult(v1,v1);
  ok &= VerifyVectorIdentity(v,v1,gVerbose,epsilon);
  if (gVerbose)
    cout << "   compare (v = pattern^2)/pattern with pattern" << endl;
  v = pattern; v1 = pattern;
  v.Sqr();
  ElementDiv(v,v1);
  ok &= VerifyVectorValue(v,pattern,gVerbose,epsilon);
  if (gVerbose)
   Compare(v1,v);

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(3,"Binary vector element-by-element operations",ok);
}

//
//------------------------------------------------------------------------
//               Verify the norm calculation
//
void vstress_norms(Int_t vsize)
{
  if (gVerbose)
    cout << "\n---> Verify norm calculations" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 10.25;

  if ( vsize % 2 == 1 )
    Fatal("vstress_norms", "size of the vector to test must be even for this test\n");

  TVectorD v(vsize);
  TVectorD v1(v);

  if (gVerbose)
    cout << "\nAssign " << pattern << " to all the elements and check norms" << endl;
  v = pattern;
  if (gVerbose)
    cout << "  1. norm should be pattern*no_elems" << endl;
  ok &= ( v.Norm1() == pattern*v.GetNoElements() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Square of the 2. norm has got to be pattern^2 * no_elems" << endl;
  ok &= ( v.Norm2Sqr() == (pattern*pattern)*v.GetNoElements() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Inf norm should be pattern itself" << endl;
  ok &= ( v.NormInf() == pattern ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Scalar product of vector by itself is the sqr(2. vector norm)" << endl;
  ok &= ( v.Norm2Sqr() == v*v ) ? kTRUE : kFALSE;

  Double_t ap_step = 1;
  Double_t ap_a0   = -pattern;
  Int_t n = v.GetNoElements();
  if (gVerbose) {
    cout << "\nAssign the arithm progression with 1. term " << ap_a0 <<
            "\nand the difference " << ap_step << endl;
  }
  {
    for (Int_t i = v.GetLwb(); i <= v.GetUpb(); i++)
      v(i) = (i-v.GetLwb())*ap_step + ap_a0;
  }
  Int_t l = TMath::Min(TMath::Max((int)TMath::Ceil(-ap_a0/ap_step),0),n);
  Double_t norm = (2*ap_a0+(l+n-1)*ap_step)/2*(n-l) +
                  (-2*ap_a0-(l-1)*ap_step)/2*l;
  if (gVerbose)
    cout << "  1. norm should be " << norm << endl;
  ok &= ( v.Norm1() == norm ) ? kTRUE : kFALSE;
  norm = n*( (ap_a0*ap_a0)+ap_a0*ap_step*(n-1)+(ap_step*ap_step)*(n-1)*(2*n-1)/6);
  if (gVerbose) {
    cout << "  Square of the 2. norm has got to be "
            "n*[ a0^2 + a0*q*(n-1) + q^2/6*(n-1)*(2n-1) ], or " << norm << endl;
  }
  ok &= ( TMath::Abs( (v.Norm2Sqr()-norm)/norm ) < EPSILON ) ? kTRUE : kFALSE;

  norm = TMath::Max(TMath::Abs(v(v.GetLwb())),TMath::Abs(v(v.GetUpb())));
  if (gVerbose)
    cout << "  Inf norm should be max(abs(a0),abs(a0+(n-1)*q)), ie " << norm << endl;
  ok &= ( v.NormInf() == norm ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "  Scalar product of vector by itself is the sqr(2. vector norm)" << endl;
  ok &= ( v.Norm2Sqr() == v*v ) ? kTRUE : kFALSE;

#if 0
  v1.Zero();
  Compare(v,v1);  // they are not equal (of course)
#endif

  if (gVerbose)
    cout << "\nConstruct v1 to be orthogonal to v as v(n), -v(n-1), v(n-2)..." << endl;
  {
    for (Int_t i = 0; i < v1.GetNoElements(); i++)
      v1(i+v1.GetLwb()) = v(v.GetUpb()-i) * ( i % 2 == 1 ? -1 : 1 );
  }
  if (gVerbose)
    cout << "||v1|| has got to be equal ||v|| regardless of the norm def" << endl;
  ok &= ( v1.Norm1()    == v.Norm1() ) ? kTRUE : kFALSE;
  ok &= ( v1.Norm2Sqr() == v.Norm2Sqr() ) ? kTRUE : kFALSE;
  ok &= ( v1.NormInf()  == v.NormInf() ) ? kTRUE : kFALSE;
  if (gVerbose)
    cout << "But the scalar product has to be zero" << endl;
  ok &= ( v1 * v == 0 ) ? kTRUE : kFALSE;

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(4,"Vector Norms",ok);
}

//
//------------------------------------------------------------------------
//           Test operations with vectors and matrix slices
//
void vstress_matrix_slices(Int_t vsize)
{
  if (gVerbose)
    cout << "\n---> Test operations with vectors and matrix slices" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = 8.625;

  TVectorD vc(0,vsize);
  TVectorD vr(0,vsize+1);
  TMatrixD m(0,vsize,0,vsize+1);

  Int_t i,j;
  if (gVerbose)
    cout << "\nCheck modifying the matrix column-by-column" << endl;
  m = pattern;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetColLwb(); i <= m.GetColUpb(); i++) {
    TMatrixDColumn(m,i) = pattern-1;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    TMatrixDColumn(m,i) *= 2;
    vc = TMatrixDColumn(m,i);
    ok &= VerifyVectorValue(vc,2*(pattern-1),gVerbose);
    vc = TMatrixDColumn(m,i+1 > m.GetColUpb() ? m.GetColLwb() : i+1);
    ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
    TMatrixDColumn(m,i) *= 0.5;
    TMatrixDColumn(m,i) += 1;
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetColLwb(); i <= m.GetColUpb(); i++) {
    vc = pattern+1;
    TMatrixDColumn(m,i) = vc;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    {
      TMatrixDColumn mc(m,i);
      for (j = m.GetRowLwb(); j <= m.GetRowUpb(); j++)
        mc(j) *= 4;
    }
    vc = TMatrixDColumn(m,i);
    ok &= VerifyVectorValue(vc,4*(pattern+1),gVerbose,EPSILON);
    TMatrixDColumn(m,i) *= 0.25;
    TMatrixDColumn(m,i) += -1;
    vc = TMatrixDColumn(m,i);
    ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nCheck modifying the matrix row-by-row" << endl;
  m = pattern;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    TMatrixDRow(m,i) = pattern+2;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    vr = TMatrixDRow(m,i);
    ok &= VerifyVectorValue(vr,pattern+2,gVerbose,EPSILON);
    vr = TMatrixDRow(m,i+1 > m.GetRowUpb() ? m.GetRowLwb() : i+1);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    TMatrixDRow(m,i) += -2;
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  for (i = m.GetRowLwb(); i <= m.GetRowUpb(); i++) {
    vr = pattern-2;
    TMatrixDRow(m,i) = vr;
    ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
    {
      TMatrixDRow mr(m,i);
      for (j = m.GetColLwb(); j <= m.GetColUpb(); j++)
        mr(j) *= 8;
    }
    vr = TMatrixDRow(m,i);
    ok &= VerifyVectorValue(vr,8*(pattern-2),gVerbose,EPSILON);
    TMatrixDRow(m,i) *= 1./8;
    TMatrixDRow(m,i) += 2;
    vr = TMatrixDRow(m,i);
    ok &= VerifyVectorValue(vr,pattern,gVerbose,EPSILON);
    ok &= ( m == pattern ) ? kTRUE : kFALSE;
  }

  if (gVerbose)
    cout << "\nCheck modifying the matrix diagonal" << endl;
  m = pattern;
  TMatrixDDiag td = m;
  td = pattern-3;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  vc = TMatrixDDiag(m);
  ok &= VerifyVectorValue(vc,pattern-3,gVerbose,EPSILON);
  td += 3;
  ok &= ( m == pattern ) ? kTRUE : kFALSE;
  vc = pattern+3;
  td = vc;
  ok &= ( !( m == pattern ) && !( m != pattern ) ) ? kTRUE : kFALSE;
  {
    TMatrixDDiag md(m);
    for (j = 0; j < md.GetNdiags(); j++)
      md(j) /= 1.5;
  }
  vc = TMatrixDDiag(m);
  ok &= VerifyVectorValue(vc,(pattern+3)/1.5,gVerbose,EPSILON);
  TMatrixDDiag(m) *= 1.5;
  TMatrixDDiag(m) += -3;
  vc = TMatrixDDiag(m);
  ok &= VerifyVectorValue(vc,pattern,gVerbose,EPSILON);
  ok &= ( m == pattern ) ? kTRUE : kFALSE;

  if (gVerbose) {
    cout << "\nCheck out to see that multiplying by diagonal is column-wise"
            "\nmatrix multiplication" << endl;
  }
  TMatrixD mm(m);
  TMatrixD m1(m.GetRowLwb(),TMath::Max(m.GetRowUpb(),m.GetColUpb()),
              m.GetColLwb(),TMath::Max(m.GetRowUpb(),m.GetColUpb()));
  TVectorD vc1(vc),vc2(vc);
  for (i = m.GetRowLwb(); i < m.GetRowUpb(); i++)
    TMatrixDRow(m,i) = pattern+i;      // Make a multiplicand
  mm = m;                          // Save it

  m1 = pattern+10;
  for (i = vr.GetLwb(); i <= vr.GetUpb(); i++)
    vr(i) = i+2;
  TMatrixDDiag td2 = m1;
  td2 = vr;
  ok &= ( !(m1 == pattern+10) ) ? kTRUE : kFALSE;

  m *= TMatrixDDiag(m1);
  for (i = m.GetColLwb(); i <= m.GetColUpb(); i++) {
    vc1 = TMatrixDColumn(mm,i);
    vc1 *= vr(i);                    // Do a column-wise multiplication
    vc2 = TMatrixDColumn(m,i);
    ok &= VerifyVectorIdentity(vc1,vc2,gVerbose,EPSILON);
  }

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(5,"Matrix Slices to Vectors",ok);
}

//
//------------------------------------------------------------------------
//           Test vector I/O
//
void vstress_vector_io()
{
  if (gVerbose)
    cout << "\n---> Test vector I/O" << endl;

  Bool_t ok = kTRUE;
  const Double_t pattern = TMath::Pi();

  TFile *f = new TFile("stress-vvector.root","RECREATE");

  Char_t name[80];
  Int_t iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TVectorD v(msize);
    v = pattern;

    Double_t *pattern_array = new Double_t[msize];
    for (Int_t i = 0; i < msize; i++)
      pattern_array[i] = pattern;
    TVectorD va;
    va.Use(msize,pattern_array);

    if (verbose)
      cout << "\nWrite vector v to database" << endl;

    snprintf(name,80,"v_%d",msize);
    v.Write(name);
    snprintf(name,80,"va_%d",msize);
    va.Write(name);

    delete [] pattern_array;

    iloop--;
  }

  if (gVerbose)
    cout << "\nClose database" << endl;
  delete f;

  if (gVerbose)
    cout << "\nOpen database in read-only mode and read vector" << endl;
  TFile *f1 = new TFile("stress-vvector.root");

  iloop = gNrLoop;
  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && iloop==gNrLoop);

    TVectorD v(msize);
    v = pattern;

    snprintf(name,80,"v_%d",msize);
    TVectorD *vr  = (TVectorD*) f1->Get(name);
    snprintf(name,80,"va_%d",msize);
    TVectorD *var = (TVectorD*) f1->Get(name);

    if (verbose)
      cout << "\nRead vector should be same as original still in memory" << endl;
    ok &= ((*vr) == v)  ? kTRUE : kFALSE;
    ok &= ((*var) == v) ? kTRUE : kFALSE;

    iloop--;
  }

  delete f1;

  if (gVerbose)
    cout << "\nDone\n" << endl;
  StatusPrint(6,"Vector Persistence",ok);
}

Bool_t test_svd_expansion(const TMatrixD &A)
{
  if (gVerbose)
    cout << "\n\nSVD-decompose matrix A and check if we can compose it back\n" <<endl;

  Bool_t ok = kTRUE;

  TDecompSVD svd(A);
  if (gVerbose) {
    cout << "left factor U" <<endl;
    svd.GetU().Print();
    cout << "Vector of Singular values" <<endl;
    svd.GetSig().Print();
    cout << "right factor V" <<endl;
    svd.GetV().Print();
  }

  {
    if (gVerbose)
      cout << "\tchecking that U is orthogonal indeed, i.e., U'U=E and UU'=E" <<endl;
    const Int_t nRows = svd.GetU().GetNrows();
    const Int_t nCols = svd.GetU().GetNcols();
    TMatrixD E1(nRows,nRows); E1.UnitMatrix();
    TMatrixD E2(nCols,nCols); E2.UnitMatrix();
    TMatrixD ut(TMatrixD::kTransposed,svd.GetU());
    ok &= VerifyMatrixIdentity(ut * svd.GetU(),E2,gVerbose,100*EPSILON);
    ok &= VerifyMatrixIdentity(svd.GetU() * ut,E1,gVerbose,100*EPSILON);
  }

  {
    if (gVerbose)
      cout << "\tchecking that V is orthogonal indeed, i.e., V'V=E and VV'=E" <<endl;
    const Int_t nRows = svd.GetV().GetNrows();
    const Int_t nCols = svd.GetV().GetNcols();
    TMatrixD E1(nRows,nRows); E1.UnitMatrix();
    TMatrixD E2(nCols,nCols); E2.UnitMatrix();
    TMatrixD vt(TMatrixD::kTransposed,svd.GetV());
    ok &= VerifyMatrixIdentity(vt * svd.GetV(),E2,gVerbose,100*EPSILON);
    ok &= VerifyMatrixIdentity(svd.GetV() * vt,E1,gVerbose,100*EPSILON);
  }

  {
    if (gVerbose)
      cout << "\tchecking that U*Sig*V' is indeed A" <<endl;
    const Int_t nRows = svd.GetU().GetNrows();
    const Int_t nCols = svd.GetV().GetNcols();
    TMatrixD s(nRows,nCols);
    TMatrixDDiag diag(s); diag = svd.GetSig();
    TMatrixD vt(TMatrixD::kTransposed,svd.GetV());
    TMatrixD tmp = s * vt;
    ok &= VerifyMatrixIdentity(A,svd.GetU() * tmp,gVerbose,100*EPSILON);
    if (gVerbose) {
      cout << "U*Sig*V'" <<endl;
      (svd.GetU()*tmp).Print();
    }
  }

  return ok;
}

#ifndef __CINT__
// Make a matrix from an array (read it row-by-row)
class MakeMatrix : public TMatrixDLazy {
  const Double_t *array;
        Int_t     no_elems;
  void FillIn(TMatrixD& m) const {
    R__ASSERT( m.GetNrows() * m.GetNcols() == no_elems );
    const Double_t *ap = array;
          Double_t *mp = m.GetMatrixArray();
    for (Int_t i = 0; i < no_elems; i++)
      *mp++ = *ap++;
  }

public:
  MakeMatrix(Int_t nrows,Int_t ncols,
             const Double_t *_array,Int_t _no_elems)
    :TMatrixDLazy(nrows,ncols), array(_array), no_elems(_no_elems) {}
  MakeMatrix(Int_t row_lwb,Int_t row_upb,Int_t col_lwb,Int_t col_upb,
             const Double_t *_array,Int_t _no_elems)
    : TMatrixDLazy(row_lwb,row_upb,col_lwb,col_upb),
      array(_array), no_elems(_no_elems) {}
};
#else
TMatrixD MakeMatrix(Int_t nrows,Int_t ncols,
                    const Double_t *_array,Int_t _no_elems)
{
  TMatrixD m(nrows,ncols,_array);
  return m;
}
#endif

void astress_decomp()
{
  Bool_t ok = kTRUE;

  {
    if (gVerbose)
      cout << "\nBrandt example page 503\n" <<endl;

    Double_t array0[] = { -2,1,0,0, 2,1,0,0, 0,0,0,0, 0,0,0,0 };
    ok &= test_svd_expansion(MakeMatrix(4,4,array0,sizeof(array0)/sizeof(array0[0])));
  }

  {
    if (gVerbose)
      cout << "\nRotated by PI/2 Matrix Diag(1,4,9)\n" <<endl;

    Double_t array1[] = {0,-4,0,  1,0,0,  0,0,9 };
    ok &= test_svd_expansion(MakeMatrix(3,3,array1,sizeof(array1)/sizeof(array1[0])));
  }

  {
    if (gVerbose)
      cout << "\nExample from the Forsythe, Malcolm, Moler's book\n" <<endl;

    Double_t array2[] =
         { 1,6,11, 2,7,12, 3,8,13, 4,9,14, 5,10,15};
    ok &= test_svd_expansion(MakeMatrix(5,3,array2,sizeof(array2)/sizeof(array2[0])));
  }

  {
    if (gVerbose)
      cout << "\nExample from the Wilkinson, Reinsch's book\n" <<
              "Singular numbers are 0, 19.5959, 20, 0, 35.3270\n" <<endl;

    Double_t array3[] =
        { 22, 10,  2,   3,  7,    14,  7, 10,  0,  8,
          -1, 13, -1, -11,  3,    -3, -2, 13, -2,  4,
           9,  8,  1,  -2,  4,     9,  1, -7,  5, -1,
           2, -6,  6,   5,  1,     4,  5,  0, -2,  2 };
    ok &= test_svd_expansion(MakeMatrix(8,5,array3,sizeof(array3)/sizeof(array3[0])));
  }

  {
    if (gVerbose)
      cout << "\nExample from the Wilkinson, Reinsch's book\n" <<
              "Ordered singular numbers are Sig[21-k] = sqrt(k*(k-1))\n" <<endl;
    TMatrixD A(21,20);
    for (Int_t irow = A.GetRowLwb(); irow <= A.GetRowUpb(); irow++)
      for (Int_t icol = A.GetColLwb(); icol <= A.GetColUpb(); icol++)
        A(irow,icol) = (irow>icol ? 0 : ((irow==icol) ? 20-icol : -1));

    ok &= test_svd_expansion(A);
  }

  if (0)
  {
    if (gVerbose) {
      cout << "\nTest by W. Meier <wmeier@manu.com> to catch an obscure "
           << "bug in QR\n" <<endl;
      cout << "expect singular values to be\n"
           << "1.4666e-024   1.828427   3.828427   4.366725  7.932951\n" <<endl;
    }

    Double_t array4[] =
        {  1,  2,  0,  0,  0,
           0,  2,  3,  0,  0,
           0,  0,  0,  4,  0,
           0,  0,  0,  4,  5,
           0,  0,  0,  0,  5 };
    ok &= test_svd_expansion(MakeMatrix(5,5,array4,sizeof(array4)/sizeof(array4[0])));
  }

  {
    const TMatrixD m1 = THilbertMatrixD(5,5);
    TDecompLU lu(m1);
    ok &= VerifyMatrixIdentity(lu.GetMatrix(),m1,gVerbose,100*EPSILON);
  }

  {
    const TMatrixD m2 = THilbertMatrixD(5,5);
    const TMatrixDSym mtm(TMatrixDSym::kAtA,m2);
    TDecompChol chol(mtm);
    ok &= VerifyMatrixIdentity(chol.GetMatrix(),mtm,gVerbose,100*EPSILON);
  }

  if (gVerbose)
    cout << "\nDone" <<endl;

  StatusPrint(1,"Decomposition / Reconstruction",ok);
}

void astress_lineqn()
{
  if (gVerbose)
    cout << "\nSolve Ax=b where A is a Hilbert matrix and b(i) = sum_j Aij\n" <<endl;

  Bool_t ok = kTRUE;

  Int_t iloop = gNrLoop;
  Int_t nr    = 0;

  while (iloop >= 0) {
    const Int_t msize = gSizeA[iloop];

    const Int_t verbose = (gVerbose && nr==0 && iloop==gNrLoop);

    if (verbose)
      cout << "\nSolve Ax=b for size = " << msize <<endl;

    // Since The Hilbert matrix is accuracy "challenged", I will use a diagonaly
    // dominant one fore sizes > 100, otherwise the verification might fail

    TMatrixDSym m1 = THilbertMatrixDSym(-1,msize-2);
    TMatrixDDiag diag1(m1);
    diag1 += 1.;

    TVectorD rowsum1(-1,msize-2); rowsum1.Zero();
    TVectorD colsum1(-1,msize-2); colsum1.Zero();
    for (Int_t irow = m1.GetRowLwb(); irow <= m1.GetColUpb(); irow++) {
      for (Int_t icol = m1.GetColLwb(); icol <= m1.GetColUpb(); icol++) {
        rowsum1(irow) += m1(irow,icol);
        colsum1(icol) += m1(irow,icol);
      }
    }

    TMatrixDSym m2 = THilbertMatrixDSym(msize);
    TMatrixDDiag diag2(m2);
    diag2 += 1.;

    TVectorD rowsum2(msize); rowsum2.Zero();
    TVectorD colsum2(msize); colsum2.Zero();
    for (Int_t irow = m2.GetRowLwb(); irow <= m2.GetColUpb(); irow++) {
      for (Int_t icol = m2.GetColLwb(); icol <= m2.GetColUpb(); icol++) {
        rowsum2(irow) += m2(irow,icol);
        colsum2(icol) += m2(irow,icol);
      }
    }

    TVectorD b1(-1,msize-2);
    TVectorD b2(msize);
    {
      TDecompLU lu(m1,1.0e-20);
      b1 = rowsum1;
      lu.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
      b1 = colsum1;
      lu.TransSolve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompChol chol(m1,1.0e-20);
      b1 = rowsum1;
      chol.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
      b1 = colsum1;
      chol.TransSolve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompQRH qrh1(m1,1.0e-20);
      b1 = rowsum1;
      qrh1.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);

      TDecompQRH qrh2(m2,1.0e-20);
      b2 = colsum2;
      qrh2.TransSolve(b2);
      if (msize < 10)
        ok &= VerifyVectorValue(b2,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompSVD svd1(m1);
      b1 = rowsum1;
      svd1.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);

      TDecompSVD svd2(m2);
      b2 = colsum2;
      svd2.TransSolve(b2);
      if (msize < 10)
        ok &= VerifyVectorValue(b2,1.0,verbose,msize*EPSILON);
    }

    {
      TDecompBK bk(m1,1.0e-20);
      b1 = rowsum1;
      bk.Solve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
      b1 = colsum1;
      bk.TransSolve(b1);
      if (msize < 10)
        ok &= VerifyVectorValue(b1,1.0,verbose,msize*EPSILON);
    }

#ifndef __CINT__
    if (nr >= Int_t(1.0e+5/msize/msize)) {
#else
    if (nr >= Int_t(1.0e+3/msize/msize)) {
#endif
      nr = 0;
      iloop--;
    } else
      nr++;

    if (!ok) {
      if (gVerbose)
        cout << "Linear Equations failed for size " << msize << endl;
      break;
    }
  }

  if (gVerbose)
    cout << "\nDone" <<endl;

  StatusPrint(2,"Linear Equations",ok);
}

void astress_pseudo()
{
// The pseudo-inverse of A is found by "inverting" the SVD of A.
// To be more precise, we use SVD to solve the equation
// AX=B where B is a unit matrix.
//
// After we compute the pseudo-inverse, we verify the Moore-Penrose
// conditions: Matrix X is a pseudo-inverse of A iff
//      AXA = A
//      XAX = X
//      XA = (XA)' (i.e., XA is symmetric)
//      AX = (AX)' (i.e., AX is symmetric)

  Bool_t ok = kTRUE;

  // Allocate and fill matrix A
  enum {nrows = 4, ncols = 3};
#ifndef __CINT__
  const Double_t A_data [nrows*ncols] =
#else
  const Double_t A_data [12] =
#endif
   {0, 0, 0,
    0, 0, 0,
    1, 1, 0,
    4, 3, 0};
  TMatrixD A(nrows,ncols,A_data);

  // Perform the SVD decomposition of the transposed matrix.

  TDecompSVD svd(A);
  if (gVerbose)
    cout << "\ncondition number " << svd.Condition() << endl;

  // Compute the inverse as a solution of A*Ainv = E
  TMatrixD Ainv(nrows,nrows); Ainv.UnitMatrix();
  svd.MultiSolve(Ainv);
  Ainv.ResizeTo(ncols,nrows);
  TMatrixD Ainv2(nrows,nrows);
  svd.Invert(Ainv2);
  Ainv2.ResizeTo(ncols,nrows);
  ok &= VerifyMatrixIdentity(Ainv,Ainv2,gVerbose,EPSILON);

  if (gVerbose) {
    cout << "\nChecking the Moore-Penrose conditions for the Ainv" << endl;
    cout << "\nVerify that Ainv * A * Ainv is Ainv" << endl;
  }
  ok &= VerifyMatrixIdentity(Ainv * A * Ainv, Ainv,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nVerify that A * Ainv * A is A" << endl;
  ok &= VerifyMatrixIdentity(A * Ainv * A, A,gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nVerify that Ainv * A is symmetric" << endl;
  ok &= VerifyMatrixIdentity(Ainv * A, (Ainv*A).T(),gVerbose,EPSILON);

  if (gVerbose)
    cout << "\nVerify that A * Ainv is symmetric" << endl;
  ok &= VerifyMatrixIdentity(A * Ainv, (A*Ainv).T(),gVerbose,EPSILON);

  StatusPrint(3,"Pseudo-Inverse, Moore-Penrose",ok);
}

void astress_eigen(Int_t msize)
{
  Bool_t ok = kTRUE;

  // Check :
  // 1.  the eigen values of M' * M are the same as the singular values
  //     squared of the SVD of M .
  // 2.  M' * M  x_i =  lambda_i x_i where x_i and lambda_i are the ith
  //     eigen - vector/value .

  const TMatrixD m = THilbertMatrixD(msize,msize);

  TDecompSVD svd(m);
  TVectorD sig = svd.GetSig(); sig.Sqr();

  // Symmetric matrix EigenVector algorithm
  TMatrixDSym mtm1(TMatrixDSym::kAtA,m);
  const TMatrixDSymEigen eigen1(mtm1);
  const TVectorD eigenVal1 = eigen1.GetEigenValues();

  ok &= VerifyVectorIdentity(sig,eigenVal1,gVerbose,EPSILON);

  // Check that the eigen vectors by comparing M'M x to lambda x
  const TMatrixD mtm1x = mtm1 * eigen1.GetEigenVectors();
  TMatrixD lam_x1 = eigen1.GetEigenVectors();
  lam_x1.NormByRow(eigen1.GetEigenValues(),"");

  ok &= VerifyMatrixIdentity(mtm1x,lam_x1,gVerbose,EPSILON);

  // General matrix EigenVector algorithm tested on symmetric matrix
  TMatrixD mtm2(m,TMatrixD::kTransposeMult,m);
  const TMatrixDEigen eigen2(mtm2);
  const TVectorD eigenVal2 = eigen2.GetEigenValuesRe();

  ok &= VerifyVectorIdentity(sig,eigenVal2,gVerbose,EPSILON);

  // Check that the eigen vectors by comparing M'M x to lambda x
  const TMatrixD mtm2x = mtm2 * eigen2.GetEigenVectors();
  TMatrixD lam_x2 = eigen2.GetEigenVectors();
  lam_x2.NormByRow(eigen2.GetEigenValuesRe(),"");

  ok &= VerifyMatrixIdentity(mtm2x,lam_x2,gVerbose,EPSILON);

  // The Imaginary part of the eigenvalues should be zero
  const TVectorD eigenValIm = eigen2.GetEigenValuesIm();
  TVectorD epsVec(msize); epsVec = EPSILON;
  ok &= VerifyVectorIdentity(epsVec,eigenValIm,gVerbose,EPSILON);

  StatusPrint(4,"Eigen - Values/Vectors",ok);
}

void astress_decomp_io(Int_t msize)
{
//
//------------------------------------------------------------------------
//           Test decomposition I/O
//
  if (gVerbose)
    cout << "\n---> Test decomp I/O" << endl;

  Bool_t ok = kTRUE;

  const TMatrixDSym m = THilbertMatrixDSym(msize);
  TVectorD rowsum(msize); rowsum.Zero();
  TVectorD colsum(msize); colsum.Zero();

  for (Int_t irow = 0; irow < m.GetNrows(); irow++) {
    for (Int_t icol = 0; icol < m.GetNcols(); icol++) {
      rowsum(irow) += m(irow,icol);
      colsum(icol) += m(irow,icol);
    }
  }

  if (gVerbose)
    cout << "\nWrite decomp m to database" << endl;

  TFile *f = new TFile("stress-vdecomp.root", "RECREATE");

  TDecompLU   lu(m,1.0e-20);
  TDecompQRH  qrh(m,1.0e-20);
  TDecompChol chol(m,1.0e-20);
  TDecompSVD  svd(m);
  TDecompBK   bk(m,1.0e-20);

  lu.Write("lu");
  qrh.Write("qrh");
  chol.Write("chol");
  svd.Write("svd");
  bk.Write("bk");

  if (gVerbose)
    cout << "\nClose database" << endl;
  delete f;

  if (gVerbose)
    cout << "\nOpen database in read-only mode and read matrix" << endl;
  TFile *f1 = new TFile("stress-vdecomp.root");

  if (gVerbose)
    cout << "\nRead decompositions should create same solutions" << endl;
  {
    TDecompLU *rlu = (TDecompLU*)  f1->Get("lu");
    TVectorD b1(rowsum);
    lu.Solve(b1);
    TVectorD b2(rowsum);
    rlu->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    lu.TransSolve(b1);
    b2 = colsum;
    rlu->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompChol *rchol = (TDecompChol*) f1->Get("chol");
    TVectorD b1(rowsum);
    chol.Solve(b1);
    TVectorD b2(rowsum);
    rchol->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    chol.TransSolve(b1);
    b2 = colsum;
    rchol->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompQRH *rqrh = (TDecompQRH*) f1->Get("qrh");
    TVectorD b1(rowsum);
    qrh.Solve(b1);
    TVectorD b2(rowsum);
    rqrh->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    qrh.TransSolve(b1);
    b2 = colsum;
    rqrh->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompSVD *rsvd = (TDecompSVD*) f1->Get("svd");
    TVectorD b1(rowsum);
    svd.Solve(b1);
    TVectorD b2(rowsum);
    rsvd->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    svd.TransSolve(b1);
    b2 = colsum;
    rsvd->TransSolve(b2);
    ok &= (b1 == b2);
  }

  {
    TDecompBK *rbk = (TDecompBK*) f1->Get("bk");
    TVectorD b1(rowsum);
    bk.Solve(b1);
    TVectorD b2(rowsum);
    rbk->Solve(b2);
    ok &= (b1 == b2);
    b1 = colsum;
    bk.TransSolve(b1);
    b2 = colsum;
    rbk->TransSolve(b2);
    ok &= (b1 == b2);
  }

  delete f1;

  if (gVerbose)
    cout << "\nDone\n" << endl;

  StatusPrint(5,"Decomposition Persistence",ok);
}

void stress_backward_io()
{
  TFile *f = TFile::Open(Form("%s/linearIO.root",gInputFiles.c_str()));

  TMatrixF mf1 = THilbertMatrixF(-5,5,-5,5);
  mf1[1][2] = TMath::Pi();
  TMatrixFSym mf2 = THilbertMatrixFSym(-5,5);
  TVectorF vf_row(mf1.GetRowLwb(),mf1.GetRowUpb()); vf_row = TMatrixFRow(mf1,3);

  TMatrixF    *mf1_r    = (TMatrixF*)    f->Get("mf1");
  TMatrixFSym *mf2_r    = (TMatrixFSym*) f->Get("mf2");
  TVectorF    *vf_row_r = (TVectorF*)    f->Get("vf_row");

  Bool_t ok = kTRUE;

  ok &= ((*mf1_r) == mf1) ? kTRUE : kFALSE;
  ok &= ((*mf2_r) == mf2) ? kTRUE : kFALSE;
  ok &= ((*vf_row_r) == vf_row) ? kTRUE : kFALSE;

  TMatrixD md1 = THilbertMatrixD(-5,5,-5,5);
  md1[1][2] = TMath::Pi();
  TMatrixDSym md2 = THilbertMatrixDSym(-5,5);
  TMatrixDSparse md3 = md1;
  TVectorD vd_row(md1.GetRowLwb(),md1.GetRowUpb()); vd_row = TMatrixDRow(md1,3);

  TMatrixD       *md1_r    = (TMatrixD*)       f->Get("md1");
  TMatrixDSym    *md2_r    = (TMatrixDSym*)    f->Get("md2");
  TMatrixDSparse *md3_r    = (TMatrixDSparse*) f->Get("md3");
  TVectorD       *vd_row_r = (TVectorD*)       f->Get("vd_row");

  ok &= ((*md1_r) == md1) ? kTRUE : kFALSE;
  ok &= ((*md2_r) == md2) ? kTRUE : kFALSE;
  ok &= ((*md3_r) == md3) ? kTRUE : kFALSE;
  ok &= ((*vd_row_r) == vd_row) ? kTRUE : kFALSE;

  StatusPrint(1,"Streamers",ok);
}



// data structure for one point
typedef struct {
   Double_t x,y,z,theta,phi;      // Initial track position and direction
   Int_t    nbound;               // Number of boundaries crossed until exit
   Float_t  length;               // Total length up to exit
   Float_t  safe;                 // Safety distance for the initial location
   Float_t  rad;                  // Number of radiation lengths up to exit
} p_t;
p_t p;
  
const Int_t NG = 33;
const char *exps[NG] = {"aleph",  
                        "barres",
                        "felix",
                        "phenix",
                        "chambers",
                        "p326",
                        "bes",
                        "dubna",
                        "ganil",
                        "e907",
                        "phobos2",
                        "hermes",
                        "na35",
                        "na47",
                        "na49",
                        "wa91",
                        "sdc",
                        "integral",
                        "ams", 
                        "brahms",
                        "gem",
                        "tesla",
                        "btev",
                        "cdf",  
                        "hades2", 
                        "lhcbfull",
                        "star", 
                        "sld",   
                        "cms",   
                        "alice3",
                        "babar2", 
                        "belle",
                        "atlas" 
};
const Int_t versions[NG] =  {4, //aleph
                             3, //barres
                             3, //felix
                             3, //phenix
                             3, //chambers
                             3, //p326
                             3, //bes
                             3, //dubna
                             3, //ganil
                             3, //e907
                             4, //phobos2
                             3, //hermes
                             3, //na35
                             3, //na47
                             3, //na49
                             3, //wa91
                             3, //sdc
                             3, //integral
                             3, //ams
                             3, //brahms
                             3, //gem
                             3, //tesla
                             3, //btev
                             4, //cdf
                             4, //hades2
                             3, //lhcbfull
                             3, //star
                             3, //sld
                             3, //cms
                             4, //alice3
                             3, //babar2
                             3, //belle
                             4}; //atlas

#if VERBOSE
// The timings below are on my machine PIV 3GHz
const Double_t cp_brun[NG] = {1.9,  //aleph
                              0.1,  //barres
                              0.12, //felix
                              0.62, //phenix
                              0.1,  //chambers
                              0.19, //p326
                              1.2,  //bes
                              0.12, //dubna
                              0.11, //ganil
                              0.17, //e907
                              0.22, //phobos2
                              0.24, //hermes
                              0.14, //na35
                              0.21, //na47
                              0.23, //na49
                              0.16, //wa91
                              0.17, //sdc
                              0.63, //integral
                              0.9,  //ams
                              1.1,  //brahms
                              1.8,  //gem
                              1.5,  //tesla
                              1.6,  //btev
                              2.2,  //cdf
                              1.2,  //hades2
                              1.6,  //lhcbfull
                              2.7,  //star
                              3.3,  //sld
                              7.5,  //cms
                              8.0,  //alice2
                             19.6,  //babar2
                             24.1,  //belle
                             26.7}; //atlas

#endif

// Bounding boxes for experiments
Double_t boxes[NG][3] = {{600,600,500},     // aleph
                         {100,100,220},     // barres
                         {200,200,12000},   // felix
                         {750,750,1000},    // phenix
                         {500,500,500},     // chambers
                         {201,201,26000},   // p326
                         {400,400,240},     // bes
                         {500,500,2000},    // dubna
                         {500,500,500},     // ganil
                         {250,250,2000},    // e907
                         {400,40,520},      // phobos2
                         {250,250,770},     // hermes
                         {310,160,1500},    // na35
                         {750,500,3000},    // na47
                         {600,200,2000},    // na49
                         {175,325,680},     // wa91
                         {1400,1400,2100},  // sdc
                         {100,100,200},     // integral
                         {200,200,200},     // ams
                         {50,50,50},        // brahms
                         {2000,2000,5000},  // gem
                         {1500,1500,1500},  // tesla
                         {600,475,1270},    // btev
                         {500,500,500},     // cdf
                         {250,250,200},     // hades2
                         {6700,5000,19000}, // lhcbfull
                         {350,350,350},     // star
                         {500,500,500},     // sld
                         {800,800,1000},    // cms
                         {400,400,400},     // alice2
                         {300,300,400},     // babar2
                         {440,440,538},     // belle
                         {1000,1000,1500}   // atlas
};                     

// Total and reference times
Double_t tpstot = 0;
Double_t tpsref = 112.1; //time including the generation of the ref files
Bool_t testfailed = kFALSE;
                         
Int_t iexp[NG];
Bool_t gen_ref=kFALSE;
void FindRad(Double_t x, Double_t y, Double_t z,Double_t theta, Double_t phi, Int_t &nbound, Float_t &length, Float_t &safe, Float_t &rad, Bool_t verbose=kFALSE);
void ReadRef(Int_t kexp);
void WriteRef(Int_t kexp);
void InspectRef(const char *exp="alice", Int_t vers=3);

Double_t stressGeometry(const char *exp, Bool_t generate_ref)
 {
   TGeoManager::SetVerboseLevel(0);
   gen_ref = generate_ref;
   gErrorIgnoreLevel = 10;
   
   #if VERBOSE
   fprintf(stderr,"******************************************************************\n");
   fprintf(stderr,"* STRESS GEOMETRY\n");
   #endif

   TString opt = exp;
   opt.ToLower();
   Bool_t all = kFALSE;
   if (opt.Contains("*")) all = kTRUE;
   Int_t i;
   for (i=0; i<NG; i++) {
      if (all) {
         iexp[i] = 1;
         continue;
      }
      if (opt.Contains(exps[i])) iexp[i] = 1;
      else                       iexp[i] = 0;
   }
   iexp[NG-1]=0;
   TFile::SetCacheFileDir(".");
   TString fname;
   for (i=0; i<NG; i++) {
      if (!iexp[i]) continue;
      fname = TString::Format("%s.root", exps[i]);
      if (gGeoManager) {
         delete gGeoManager;
         gGeoManager = 0;
      }   
      TGeoManager::Import(Form("%s/%s",gInputFiles.c_str(),fname.Data()));
      if (!gGeoManager) return 0.0;
         
      fname = TString::Format("%s/%s_ref_%d.root", gInputFiles.c_str(),exps[i],versions[i]);
      
      if (gen_ref || !TFile::Open(Form("%s/%s_ref_%d.root",gInputFiles.c_str(),exps[i],versions[i]))) {
         if (!gen_ref) fprintf(stderr,"File: %s does not exist, generating it\n", fname.Data());
         else               
         {
#if VERBOSE          
          fprintf(stderr,"Generating reference file %s\n", fname.Data());
#endif          
         }
         WriteRef(i);
      }
   
      ReadRef(i);
   }   
   if (all && tpstot>0) {
#if VERBOSE      
      Float_t rootmarks = 800*tpsref/tpstot;     
      fprintf(stderr,"******************************************************************\n");
      if (testfailed) fprintf(stderr,"*  stressGeometry found bad points ............. FAILED\n");
      else          fprintf(stderr,"*  stressGeometry .................................. OK\n");
      fprintf(stderr,"******************************************************************\n");
      fprintf(stderr,"*  CPU time in ReadRef = %6.2f seconds\n",tpstot);
      fprintf(stderr,"*  streeGeometry * ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),gROOT->GetVersionDate(),gROOT->GetVersionTime());
      fprintf(stderr,"******************************************************************\n");
#endif      
      return tpstot;
   }
   return 0.0;
}

void ReadRef(Int_t kexp) {
   TStopwatch sw;
   TString fname;
   TFile *f = 0;
   //use ref_[version[i]] files
   if (!gen_ref)
      fname = TString::Format("%s/%s_ref_%d.root", gInputFiles.c_str(),exps[kexp],versions[kexp]);
   else
      fname.Format("%s/%s_ref_%d.root", gInputFiles.c_str(),exps[kexp],versions[kexp]);
   
   f = TFile::Open(fname);
   if (!f) {
      fprintf(stderr,"Reference file %s not found ! Skipping.\n", fname.Data());
      return;
   }   
   // fprintf(stderr,"Reference file %s found\n", fname.Data());
   fname = TString::Format("%s_diff.root", exps[kexp]);
   TFile fdiff(fname,"RECREATE");
   TTree *TD = new TTree("TD","TGeo stress diff");
   TD->Branch("p",&p.x,"x/D:y/D:z/D:theta/D:phi/D:rad[4]/F");
   TTree *T = (TTree*)f->Get("T");
   T->SetBranchAddress("p",&p.x);
   Long64_t nentries = T->GetEntries();
   TVectorD *vref = (TVectorD *)T->GetUserInfo()->At(0);
   if (!vref) {
      fprintf(stderr," ERROR: User info not found, regenerate reference file\n");
      return;
   }   
   TVectorD vect(4);
   TVectorD vect_ref = *vref;
   Int_t nbound;
   Float_t length, safe, rad;
   Float_t diff;
   Float_t diffmax = 0.01;  // percent of rad!
   Int_t nbad = 0;
   vect(0) = 0;//gGeoManager->Weight(0.01, "va");
   for (Long64_t i=0;i<nentries;i++) {
      T->GetEntry(i);
      nbound = 0;
      length = 0.;
      safe = 0.;
      rad = 0.;
      FindRad(p.x,p.y,p.z, p.theta, p.phi, nbound, length, safe, rad);
      vect(1) += Double_t(nbound);
      vect(2) += length;
      vect(3) += rad;
      diff = 0;
      diff += TMath::Abs(length-p.length);
      diff += TMath::Abs(safe-p.safe);
      diff += TMath::Abs(rad-p.rad);
      if (((p.rad>0) && (TMath::Abs(rad-p.rad)/p.rad)>diffmax) || 
           TMath::Abs(nbound-p.nbound)>100) {
         nbad++;
         if (nbad < 10) {
            fprintf(stderr," ==>Point %lld differs with diff = %g, x=%g, y=%g, z=%g\n",i,diff,p.x,p.y,p.z);
            fprintf(stderr,"    p.nbound=%d, p.length=%g, p.safe=%g, p.rad=%g\n",
                        p.nbound,p.length,p.safe,p.rad);
            fprintf(stderr,"      nbound=%d,   length=%g,   safe=%g,   rad=%g\n",
                        nbound,length,safe,rad);
         }
         TD->Fill();
         p.nbound = nbound;
         p.length = length;
         p.safe   = safe;
         p.rad    = rad;
         TD->Fill();
      }    
   }
   diff = 0.;
   //for (Int_t j=1; j<4; j++) diff += TMath::Abs(vect_ref(j)-vect(j));
   diff += TMath::Abs(vect_ref(3)-vect(3))/vect_ref(3);
   if (diff > diffmax) {
//      fprintf(stderr,"Total weight=%g   ref=%g\n", vect(0), vect_ref(0));
      fprintf(stderr,"Total nbound=%g   ref=%g\n", vect(1), vect_ref(1));
      fprintf(stderr,"Total length=%g   ref=%g\n", vect(2), vect_ref(2));
      fprintf(stderr,"Total    rad=%g   ref=%g\n", vect(3), vect_ref(3));
      nbad++;  
   }   
      
   if (nbad) {
      testfailed = kTRUE;
      TD->AutoSave();
      TD->Print();
   }   
   delete TD;
   delete f;
   
   Double_t cp = sw.CpuTime();
   tpstot += cp;
   if (nbad > 0) 
   {
    fprintf(stderr,"*     stress %-15s  found %5d bad points ............. failed\n",exps[kexp],nbad);
  }
   else
   {
#if VERBOSE
             fprintf(stderr,"*     stress %-15s: time/ref = %6.2f/%6.2f............ OK\n",exps[kexp],cp,cp_brun[kexp]);
#endif             
           }
}

void WriteRef(Int_t kexp) {
   TRandom3 r;
//   Double_t theta, phi;
   Double_t point[3];
   TVectorD vect(4);
   TGeoShape *top = gGeoManager->GetMasterVolume()->GetShape();
//   TGeoBBox *box = (TGeoBBox*)top;
   Double_t xmax = boxes[kexp][0]; //box->GetDX(); // 300;
   Double_t ymax = boxes[kexp][1]; //box->GetDY(); // 300;
   Double_t zmax = boxes[kexp][2]; //box->GetDZ(); // 500;
   TString fname(TString::Format("%s/%s_ref_%d.root", gInputFiles.c_str(), exps[kexp], versions[kexp]));
   TFile f(fname,"recreate");
   TTree *T = new TTree("T","TGeo stress");
   T->Branch("p",&p.x,"x/D:y/D:z/D:theta/D:phi/D:nbound/I:length/F:safe/F:rad/F");
   T->GetUserInfo()->Add(&vect);
   Long64_t Npoints = 10000;
   Long64_t i = 0;
   vect(0) = 0; //gGeoManager->Weight(0.01, "va");
   while (i<Npoints) {
      p.x  = r.Uniform(-xmax,xmax);
      p.y  = r.Uniform(-ymax,ymax);
      p.z  = r.Uniform(-zmax,zmax);
      point[0] = p.x;
      point[1] = p.y;
      point[2] = p.z;
      if (top->Contains(point)) {
         p.phi   =  2*TMath::Pi()*r.Rndm();
         p.theta = TMath::ACos(1.-2.*r.Rndm());
         FindRad(p.x,p.y,p.z, p.theta, p.phi, p.nbound, p.length, p.safe, p.rad);
         vect(1) += Double_t(p.nbound);
         vect(2) += p.length;
         vect(3) += p.rad;
         T->Fill();
         i++;
      }
   }   
   T->AutoSave();
   T->GetUserInfo()->Remove(&vect);
//   T->Print();
   delete T;
}

void FindRad(Double_t x, Double_t y, Double_t z,Double_t theta, Double_t phi, Int_t &nbound, Float_t &length, Float_t &safe, Float_t &rad, Bool_t verbose) {
   Double_t xp  = TMath::Sin(theta)*TMath::Cos(phi);
   Double_t yp  = TMath::Sin(theta)*TMath::Sin(phi);
   Double_t zp  = TMath::Cos(theta);
   Double_t snext;
   TString path;
   Double_t pt[3];
   Double_t loc[3];
   Double_t epsil = 1.E-2;
   Double_t lastrad = 0.;
   Int_t ismall = 0;
   nbound = 0;
   length = 0.;
   safe   = 0.;
   rad    = 0.;
   TGeoMedium *med;
   TGeoShape *shape;
   TGeoNode *lastnode;
   gGeoManager->InitTrack(x,y,z,xp,yp,zp);
   if (verbose) {
      fprintf(stderr,"Track: (%15.10f,%15.10f,%15.10f,%15.10f,%15.10f,%15.10f)\n",
                       x,y,z,xp,yp,zp);
      path = gGeoManager->GetPath();
   }                    
   TGeoNode *nextnode = gGeoManager->GetCurrentNode();
   safe = gGeoManager->Safety();
   while (nextnode) {
      med = 0;
      if (nextnode) med = nextnode->GetVolume()->GetMedium();
      else return;      
      shape = nextnode->GetVolume()->GetShape();
      lastnode = nextnode;
      nextnode = gGeoManager->FindNextBoundaryAndStep();
      snext  = gGeoManager->GetStep();
      if (snext<1.e-8) {
         ismall++;
         if ((ismall<3) && (lastnode != nextnode)) {
            // First try to cross a very thin layer
            length += snext;
            nextnode = gGeoManager->FindNextBoundaryAndStep();
            snext  = gGeoManager->GetStep();
            if (snext<1.E-8) continue;
            // We managed to cross the layer
            ismall = 0;
         } else {  
            // Relocate point
            if (ismall > 3) {
               fprintf(stderr,"ERROR: Small steps in: %s shape=%s\n",gGeoManager->GetPath(), shape->ClassName());
               return;
            }   
            memcpy(pt,gGeoManager->GetCurrentPoint(),3*sizeof(Double_t));
            const Double_t *dir = gGeoManager->GetCurrentDirection();
            for (Int_t i=0;i<3;i++) pt[i] += epsil*dir[i];
            snext = epsil;
            length += snext;
            rad += lastrad*snext;
            gGeoManager->CdTop();
            nextnode = gGeoManager->FindNode(pt[0],pt[1],pt[2]);
            if (gGeoManager->IsOutside()) return;
            TGeoMatrix *mat = gGeoManager->GetCurrentMatrix();
            mat->MasterToLocal(pt,loc);
            if (!gGeoManager->GetCurrentVolume()->Contains(loc)) {
//            fprintf(stderr,"Woops - out\n");
               gGeoManager->CdUp();
               nextnode = gGeoManager->GetCurrentNode();
            }   
            continue;
         }   
      } else {
         ismall = 0;
      }      
      nbound++;
      length += snext;
      if (med) {
         Double_t radlen = med->GetMaterial()->GetRadLen();
         if (radlen>1.e-5 && radlen<1.e10) {
            lastrad = med->GetMaterial()->GetDensity()/radlen;
            rad += lastrad*snext;
         } else {
            lastrad = 0.;
         }      
         if (verbose) {
            fprintf(stderr," STEP #%d: %s\n",nbound, path.Data());
            fprintf(stderr,"    step=%g  length=%g  rad=%g %s\n", snext,length,
                   med->GetMaterial()->GetDensity()*snext/med->GetMaterial()->GetRadLen(),med->GetName());
            path =  gGeoManager->GetPath();
         }   
      }
   }   
}
  
void InspectDiff(const char* exp="alice",Long64_t ientry=-1) {
   Int_t nbound = 0;   
   Float_t length = 0.;
   Float_t safe   = 0.;
   Float_t rad    = 0.;
   TString fname(TString::Format("%s.root",exp));
   if (gSystem->AccessPathName(fname)) {
      TGeoManager::Import(Form("%s/%s",gInputFiles.c_str(),fname.Data()));
   } else {
      TGeoManager::Import(fname);
   }
   fname = TString::Format("%s_diff.root",exp);   
   TFile f(fname);
   if (f.IsZombie()) return;
   TTree *TD = (TTree*)f.Get("TD");
   TD->SetBranchAddress("p",&p.x);
   Long64_t nentries = TD->GetEntries();
   nentries = nentries>>1;
   if (ientry>=0 && ientry<nentries) {
      fprintf(stderr,"DIFFERENCE #%lld\n", ientry);
      TD->GetEntry(2*ientry);
      fprintf(stderr,"   NEW: nbound=%d  length=%g  safe=%g  rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      TD->GetEntry(2*ientry+1);
      fprintf(stderr,"   OLD: nbound=%d  length=%g  safe=%g  rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      FindRad(p.x,p.y,p.z, p.theta, p.phi, nbound,length,safe,rad, kTRUE);
      return;
   }   
   for (Long64_t i=0;i<nentries;i++) {
      fprintf(stderr,"DIFFERENCE #%lld\n", i);
      TD->GetEntry(2*i);
      fprintf(stderr,"   NEW: nbound=%d  length=%g  safe=%g rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      TD->GetEntry(2*i+1);
      fprintf(stderr,"   OLD: nbound=%d  length=%g  safe=%g rad=%g\n", p.nbound,p.length,p.safe,p.rad);
      FindRad(p.x,p.y,p.z, p.theta, p.phi, nbound,length,safe,rad, kTRUE);
   }
}   

void InspectRef(const char *exp, Int_t vers) {
// Inspect current reference.
   TString fname(TString::Format("%s_ref_%d.root", exp, vers));
   if (gSystem->AccessPathName(fname)) {
      fprintf(stderr,"ERROR: file %s does not exist\n", fname.Data());
      return;
   }
   TFile f(fname);
   if (f.IsZombie()) return;
   TTree *T = (TTree*)f.Get("T");
   Long64_t nentries = T->GetEntries();
   fname.Format("Stress test for %s geometry", exp);
   TCanvas *c = new TCanvas("stress", fname,700,800);
   c->Divide(2,2,0.005,0.005);
   c->cd(1);
   gPad->SetLogy();
   T->Draw("p.nbound","","", nentries, 0);
   c->cd(2);
   gPad->SetLogy();
   T->Draw("p.length","","", nentries, 0);
   c->cd(3);
   gPad->SetLogy();
   T->Draw("p.safe","","", nentries, 0);
   c->cd(4);
   gPad->SetLogy();
   T->Draw("p.rad","","", nentries, 0);
   c->cd(0);
   c->SetFillColor(kYellow);
   TVectorD *vref = (TVectorD *)T->GetUserInfo()->At(0);
   TVectorD vect = *vref;
   fprintf(stderr,"=====================================\n");
//   fprintf(stderr,"Total weight:  %g [kg]\n", vect(0));
   fprintf(stderr,"Total nbound:  %g boundaries crossed\n", vect(1));
   fprintf(stderr,"Total length:  %g [m]\n", 0.01*vect(2));
   fprintf(stderr,"Total nradlen: %f\n", vect(3));   
   fprintf(stderr,"=====================================\n");
}

Int_t npeaks;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}
Double_t fpeaks2(Double_t *x, Double_t *par) {
   Double_t result = 0.1;
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm   = par[5*p+0];
      Double_t mean1  = par[5*p+1];
      Double_t sigma1 = par[5*p+2];
      Double_t mean2  = par[5*p+3];
      Double_t sigma2 = par[5*p+4];
      result += norm*TMath::Gaus(x[0],mean1,sigma1)*TMath::Gaus(x[1],mean2,sigma2);
   }
   return result;
}
void findPeaks(Int_t pmin, Int_t pmax, Int_t &nfound, Int_t &ngood, Int_t &nghost) {
   npeaks = (Int_t)gRandom->Uniform(pmin,pmax);
   Int_t nbins = 500;
   Double_t dxbins = 2;
   TH1F *h = new TH1F("h","test",nbins,0,nbins*dxbins);
   //generate n peaks at random
   Double_t par[3000];
   par[0] = 0.8;
   par[1] = -0.6/1000;
   Int_t p,pf;
   for (p=0;p<npeaks;p++) {
      par[3*p+2] = 1;
      par[3*p+3] = 10+gRandom->Rndm()*(nbins-20)*dxbins;
      par[3*p+4] = 3+2*gRandom->Rndm();
   }
   TF1 *f = new TF1("f",fpeaks,0,nbins*dxbins,2+3*npeaks);
   f->SetNpx(1000);
   f->SetParameters(par);
   h->FillRandom("f",200000);
   TSpectrum *s = new TSpectrum(4*npeaks);
   nfound = s->Search(h,2,"goff");
   //Search found peaks
   ngood = 0;
   Float_t *xpeaks = s->GetPositionX();
   for (p=0;p<npeaks;p++) {
      for (Int_t pf=0;pf<nfound;pf++) {
         Double_t dx = TMath::Abs(xpeaks[pf] - par[3*p+3]);
         if (dx <dxbins) ngood++;
      }
   }
   //Search ghost peaks
   nghost = 0;
   for (pf=0;pf<nfound;pf++) {
      Int_t nf=0;
      for (Int_t p=0;p<npeaks;p++) {
         Double_t dx = TMath::Abs(xpeaks[pf] - par[3*p+3]);
         if (dx <dxbins) nf++;
      }
      if (nf == 0) nghost++;
   }
   delete f;
   delete h;
   delete s;
}

void stress1(Int_t ntimes) {
   Int_t pmin = 5;
   Int_t pmax = 55;
   TCanvas *c1 = new TCanvas("c1","Spectrum results",10,10,800,800);
   c1->Divide(2,2);
   gStyle->SetOptFit();
   TH1F *hpeaks = new TH1F("hpeaks","Number of peaks",pmax-pmin,pmin,pmax);
   TH1F *hfound = new TH1F("hfound","% peak founds",100,0,100);
   TH1F *hgood  = new TH1F("hgood", "% good peaks",100,0,100);
   TH1F *hghost = new TH1F("hghost","% ghost peaks",100,0,100);
   Int_t nfound,ngood,nghost;
   for (Int_t i=0;i<ntimes;i++) {
      findPeaks(pmin,pmax,nfound,ngood,nghost);
      hpeaks->Fill(npeaks);
      hfound->Fill(100*Double_t(nfound)/Double_t(npeaks));
      hgood->Fill(100*Double_t(ngood)/Double_t(npeaks));
      hghost->Fill(100*Double_t(nghost)/Double_t(npeaks));
      //printf("npeaks = %d, nfound = %d, ngood = %d, nghost = %d\n",npeaks,nfound,ngood,nghost);
   }
   c1->cd(1);
   hpeaks->Fit("pol1","lq");
   c1->cd(2);
   hfound->Fit("gaus","lq");
   c1->cd(3);
   hgood->Fit("gaus","lq");
   c1->cd(4);
   hghost->Fit("gaus","lq","",0,30);
   c1->cd();
   Double_t p1  = hfound->GetFunction("gaus")->GetParameter(1);
   Double_t ep1 = hfound->GetFunction("gaus")->GetParError(1);
   Double_t p2  = hgood->GetFunction("gaus")->GetParameter(1);
   Double_t ep2 = hgood->GetFunction("gaus")->GetParError(1);
   Double_t p3  = hghost->GetFunction("gaus")->GetParameter(1);
   Double_t ep3 = hghost->GetFunction("gaus")->GetParError(1);
   Double_t p1ref = 70.21; //ref numbers obtained with ntimes=1000
   Double_t p2ref = 65.03;
   Double_t p3ref =  8.54;
      
   //printf("p1=%g+-%g, p2=%g+-%g, p3=%g+-%g\n",p1,ep1,p2,ep2,p3,ep3);

   char sok[20];
   if (TMath::Abs(p1ref-p1) < 2*ep1 && TMath::Abs(p2ref-p2) < 2*ep2  && TMath::Abs(p3ref-p3) < 2*ep3 ) {
      snprintf(sok,20,"OK");
   } else {
      snprintf(sok,20,"failed");
   }
#if VERBOSE   
   printf("Peak1 : found =%6.2f/%6.2f, good =%6.2f/%6.2f, ghost =%5.2f/%5.2f,--- %s\n",
          p1,p1ref,p2,p2ref,p3,p3ref,sok);
#endif

}
void stress2(Int_t np2) {
   npeaks = np2;
   TRandom r;
   Int_t nbinsx = 200;
   Int_t nbinsy = 200;
   Double_t xmin   = 0;
   Double_t xmax   = (Double_t)nbinsx;
   Double_t ymin   = 0;
   Double_t ymax   = (Double_t)nbinsy;
   Double_t dx = (xmax-xmin)/nbinsx;
   Double_t dy = (ymax-ymin)/nbinsy;
   TH2F *h2 = new TH2F("h2","test",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   h2->SetStats(0);
   //generate n peaks at random
   Double_t par[3000];
   Int_t p;
   for (p=0;p<npeaks;p++) {
      par[5*p+0] = r.Uniform(0.2,1);
      par[5*p+1] = r.Uniform(xmin,xmax);
      par[5*p+2] = r.Uniform(dx,5*dx);
      par[5*p+3] = r.Uniform(ymin,ymax);
      par[5*p+4] = r.Uniform(dy,5*dy);
   }
   TF2 *f2 = new TF2("f2",fpeaks2,xmin,xmax,ymin,ymax,5*npeaks);
   f2->SetNpx(100);
   f2->SetNpy(100);
   f2->SetParameters(par);
   h2->FillRandom("f2",500000);
   //now the real stuff
   TSpectrum2 *s = new TSpectrum2(2*npeaks);
   Int_t nfound = s->Search(h2,2,"goff noMarkov");
   
   //searching good and ghost peaks (approximation)
   Int_t pf,ngood = 0;
   Float_t *xpeaks = s->GetPositionX();
   Float_t *ypeaks = s->GetPositionY();
   for (p=0;p<npeaks;p++) {
      for (Int_t pf=0;pf<nfound;pf++) {
         Double_t diffx = TMath::Abs(xpeaks[pf] - par[5*p+1]);
         Double_t diffy = TMath::Abs(ypeaks[pf] - par[5*p+3]);
         if (diffx < 2*dx && diffy < 2*dy) ngood++;
      }
   }
   if (ngood > nfound) ngood = nfound;
   //Search ghost peaks (approximation)
   Int_t nghost = 0;
   for (pf=0;pf<nfound;pf++) {
      Int_t nf=0;
      for (Int_t p=0;p<npeaks;p++) {
         Double_t diffx = TMath::Abs(xpeaks[pf] - par[5*p+1]);
         Double_t diffy = TMath::Abs(ypeaks[pf] - par[5*p+3]);
         if (diffx < 2*dx && diffy < 2*dy) nf++;
      }
      if (nf == 0) nghost++;
   }
   
   delete s;
   delete f2;
   delete h2;
   Int_t nfoundRef = 163;
   Int_t ngoodRef  = 163;
   Int_t nghostRef = 8;
   char sok[20];
   if (  TMath::Abs(nfound - nfoundRef) < 5
      && TMath::Abs(ngood - ngoodRef) < 5
      && TMath::Abs(nghost - nghostRef) < 5)  {
      snprintf(sok,20,"OK");
   } else {
      snprintf(sok,20,"failed");
   }
#if VERBOSE   
   printf("Peak2 : found =%d/%d, good =%d, ghost =%2d,---------------------------- %s\n",
          nfound,npeaks,ngood,nghost,sok);
#endif   
}
   
Double_t stressSpectrum(Int_t ntimes) 
{
#if VERBOSE  
   cout << "****************************************************************************" <<endl;
   cout << "*  Starting  stress S P E C T R U M                                        *" <<endl;
   cout << "****************************************************************************" <<endl;
#endif   
   gBenchmark->Start("stressSpectrum");
   stress1(ntimes);
   stress2(300);
   gBenchmark->Stop ("stressSpectrum");
   Double_t ct = gBenchmark->GetCpuTime("stressSpectrum");
#if VERBOSE
   Double_t reftime100 = 19.04; //pcbrun compiled
  Double_t rootmarks = 800*reftime100*ntimes/(100*ct);
   printf("****************************************************************************\n");

   gBenchmark->Print("stressSpectrum");
   printf("****************************************************************************\n");
   printf("*  stressSpectrum * ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
         gROOT->GetVersionDate(),gROOT->GetVersionTime());
   printf("****************************************************************************\n");
#endif   
   return ct;
}
   
int gPrintSubBench = 0;

Double_t ntotin=0, ntotout=0;

Double_t stressMix(Int_t nevent, Int_t style,
            Int_t printSubBenchmark, UInt_t portion)
{
   //Main control function invoking all test programs
   
   gPrintSubBench = printSubBenchmark;
   
   if (nevent < 11) nevent = 11; // must have at least 10 events
   //Delete all possible objects in memory (to execute stress several times)
   gROOT->GetListOfFunctions()->Delete();
   gROOT->GetList()->Delete();

#if VERBOSE
   printf("******************************************************************\n");
   printf("*  Starting  R O O T - S T R E S S test suite with %d events\n",nevent);
   printf("******************************************************************\n");
#endif   
   // select the branch style
   TTree::SetBranchStyle(style);

   //Run the standard test suite
   gBenchmark->Start("stress");
   if (portion&1) stress1();
   if (portion&2) stress2();
   if (portion&4) stress3();
   if (portion&8) stress4();
   if (portion&16) stress5();
   if (portion&32) stress6();
   if (portion&64) stress7();
   if (portion&128) stress8(nevent);
   if (portion&256) stress9();
   if (portion&512) stress10();
   if (portion&1024) stress11();
   if (portion&2048) stress12(12);
   if (portion&4096) stress13();
   if (portion&8192) stress14();
   if (portion&16384) stress15();
   if (portion&32768) stress16();
   gBenchmark->Stop("stress");

   Float_t ct = gBenchmark->GetCpuTime("stress");

#if VERBOSE
   printf("******************************************************************\n");
   Float_t mbtot = (Float_t)(ntotin+ntotout)/1000000.;
   Float_t mbin  = (Float_t)ntotin/1000000.;
   Float_t mbout = (Float_t)ntotout/1000000.;
   printf("stress    : Total I/O =%7.1f Mbytes, I =%7.1f, O =%6.1f\n",mbtot,mbin,mbout);
   Float_t mbin1  = (Float_t)(TFile::GetFileBytesRead()/1000000.);
   Float_t mbout1 = (Float_t)(TFile::GetFileBytesWritten()/1000000.);
   Float_t mbtot1 = mbin1+mbout1;
   printf("stress    : Compr I/O =%7.1f Mbytes, I =%7.1f, O =%6.1f\n",mbtot1,mbin1,mbout1);
   gBenchmark->Print("stress");
   Float_t cp_brun_30   = 31.03;  //The difference is essentially coming from stress16
   Float_t cp_brun_1000 = 84.30;
   Float_t cp_brun = cp_brun_1000 - (cp_brun_1000 - cp_brun_30)*(1000-nevent)/(1000-30);
   Float_t rootmarks = 600*cp_brun/ct;
   printf("******************************************************************\n");
   printf("*  stress * ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),gROOT->GetVersionDate(),gROOT->GetVersionTime());
   printf("******************************************************************\n");
#endif
   return ct;
}

//_______________________________________________________________
Double_t f1int(Double_t *x, Double_t *p)
{
   //Compute a function sum of 3 gaussians
   Double_t e1 = (x[0]-p[1])/p[2];
   Double_t e2 = (x[0]-p[4])/p[5];
   Double_t e3 = (x[0]-p[7])/p[8];
   Double_t f  = p[0]*TMath::Exp(-0.5*e1*e1)
                +p[3]*TMath::Exp(-0.5*e2*e2)
                +p[6]*TMath::Exp(-0.5*e3*e3);
   return f;
}

//_______________________________________________________________
void Bprint(Int_t id, const char *title)
{
#if VERBOSE
  // Print test program number and its title
   const Int_t kMAX = 65;
   char header[80];
   snprintf(header,80,"Test %2d : %s",id,title);
   Int_t nch = strlen(header);
   for (Int_t i=nch;i<kMAX;i++) header[i] = '.';
   header[kMAX] = 0;
   header[kMAX-1] = ' ';
   printf("%s",header);
#endif   
}

//_______________________________________________________________
void stress1()
{
   //Generate two functions supposed to produce the same result
   //One function "f1form" will be computed by the TFormula class
   //The second function "f1int" will be
   //   - compiled when running in batch mode
   //   - interpreted by CINT when running in interactive mode

   Bprint(1,"Functions, Random Numbers, Histogram Fits");

   //Start with a function inline expression (managed by TFormula)
   Double_t f1params[9] = {100,-3,3,60,0,0.5,40,4,0.7};
   TF1 *f1form = new TF1("f1form","gaus(0)+gaus(3)+gaus(6)",-10,10);
   f1form->SetParameters(f1params);

   //Create an histogram and fill it randomly with f1form
   gRandom->SetSeed(65539);
   TH1F *h1form = new TH1F("h1form","distribution from f1form",100,-10,10);
   TH1F *h1diff = (TH1F*)h1form->Clone();
   h1diff->SetName("h1diff");
   h1form->FillRandom("f1form",10000);

   //Fit h1form with original function f1form
   h1form->Fit("f1form","q0");

   //same operation with an interpreted function f1int
   TF1 *f1 = new TF1("f1int",f1int,-10,10,9);
   f1->SetParameters(f1params);

   //Create an histogram and fill it randomly with f1int
   gRandom->SetSeed(65539); //make sure we start with the same random numbers
   TH1F *h1int = new TH1F("h1int","distribution from f1int",100,-10,10);
   h1int->FillRandom("f1int",10000);

   //Fit h1int with original function f1int
   h1int->Fit("f1int","q0");

   //The difference between the two histograms must be null
   h1diff->Add(h1form, h1int, 1, -1);
   Double_t hdiff = h1diff->Integral(0,101);

   //Compare fitted parameters and value of integral of f1form in [-8,6]
   Int_t npar = f1form->GetNpar();
   Double_t pdiff, pdifftot = 0;
   for (Int_t i=0;i<npar;i++) {
      pdiff = (f1form->GetParameter(i) - f1->GetParameter(i))/f1form->GetParameter(i);
      pdifftot += TMath::Abs(pdiff);
   }
   // The integral in the range [-8,6] must be = 1923.74578
   Double_t rint = TMath::Abs(f1form->Integral(-8,6) - 1923.74578);

   //Some slight differences are authorized to take into account
   //different math libraries used by the compiler, CINT and TFormula
   Bool_t OK = kTRUE;
   if (hdiff > 0.1 || pdifftot > 2.e-3 || rint > 10) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s hdiff=%g, pdifftot=%g, rint=%g\n"," ",hdiff,pdifftot,rint);
   }
   if (gPrintSubBench) { printf("Test  1 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
   //Save all objects in a Root file (will be checked by stress2)
   TFile local("stress.root","recreate");
   f1form->Write();
   f1->Write();
   h1form->Write();
   h1int->Write();
   ntotout += local.GetBytesWritten();
   //do not close the file. should be done by the destructor automatically
   delete h1int;
   delete h1form;
   delete h1diff;
}

//_______________________________________________________________
void stress2()
{
   //check length and compression factor in stress.root
   Bprint(2,"Check size & compression factor of a Root file");
   TFile f("stress.root");
   Long64_t last = f.GetEND();
   Float_t comp = f.GetCompressionFactor();

   Bool_t OK = kTRUE;
   Long64_t lastgood = 9428;
   if (last <lastgood-200 || last > lastgood+200 || comp <2.0 || comp > 2.4) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s last =%lld, comp=%f\n"," ",last,comp);
   }
   if (gPrintSubBench) { printf("Test  2 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress3()
{
   //Open stress.root, read all objects, save 10 times and purge
   //This function tests the generation and reuse of gaps in files

   Bprint(3,"Purge, Reuse of gaps in TFile");
   TFile f("stress.root","update");
   f.ReadAll();
   for (Int_t i=0;i<10;i++) {
      f.Write();
   }
   f.Purge();
   f.Write();

   //check length and compression level in stress.root
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();
   Long64_t last = f.GetEND();
   Float_t comp = f.GetCompressionFactor();
   Bool_t OK = kTRUE;
   Long64_t lastgood = 49203;
   if (last <lastgood-900 || last > lastgood+900 || comp <1.8 || comp > 2.4) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s last =%lld, comp=%f\n"," ",last,comp);
   }
   if (gPrintSubBench) { printf("Test  3 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress4()
{
// Test of 2-d histograms, functions, 2-d fits

   Bprint(4,"Test of 2-d histograms, functions, 2-d fits");

   Double_t f2params[15] = {100,-3,3,-3,3,160,0,0.8,0,0.9,40,4,0.7,4,0.7};
   TF2 *f2form = new TF2("f2form","xygaus(0)+xygaus(5)+xygaus(10)",-10,10,-10,10);
   f2form->SetParameters(f2params);

   //Create an histogram and fill it randomly with f2form
   gRandom->SetSeed(65539);
   TH2F *h2form = new TH2F("h2form","distribution from f2form",40,-10,10,40,-10,10);
   Int_t nentries = 100000;
   h2form->FillRandom("f2form",nentries);
   //Fit h2form with original function f2form
   Float_t ratio = 4*nentries/100000;
   f2params[ 0] *= ratio;
   f2params[ 5] *= ratio;
   f2params[10] *= ratio;
   f2form->SetParameters(f2params);
   h2form->Fit("f2form","q0");
   //Update stress.root
   TFile f("stress.root","update");
   h2form->Write();
   f2form->Write();

   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();

   //Compare results of fit with expected parameters
   Bool_t OK = kTRUE;
   for (int k = 0; k < 3; ++k) { 
      for (int  l = 1; l < 5; ++l) { 
         int idx = k*5+l;
         Double_t dp0  = TMath::Abs((f2form->GetParameter(idx) -f2params[idx]));
         if (f2params[idx] != 0.) dp0 /=  f2params[idx];
         bool testok =  (dp0 < 5.e-2); 
         if (!testok) {
            printf("\nfailed:   ipar=%d delta=%g, par=%g, nom=%g",idx,dp0,f2form->GetParameter(idx),f2params[idx]);
         }
         OK &= testok;
      }
   }
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    printf("\ntest failed !\n");
   if (gPrintSubBench) { printf("Test  4 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress5()
{
// Test of Postscript.
// Make a complex picture. Verify number of lines on ps file
// Testing automatically the graphics package is a complex problem.
// The best way we have found is to generate a Postscript image
// of a complex canvas containing many objects.
// The number of lines in the ps file is compared with a reference run.
// A few lines (up to 2 or 3) of difference may be expected because
// Postscript works with floats. The date and time of the run are also
// different.
// You can also inspect visually the ps file with a ps viewer.

   Bprint(5,"Test graphics & Postscript");

   TCanvas *c1 = new TCanvas("stress-canvas","stress canvas",800,600);
   gROOT->LoadClass("TPostScript","Postscript");
   TPostScript ps("stress.ps",112);

   //Get objects generated in previous test
   TFile f("stress.root");
   TF1  *f1form = (TF1*)f.Get("f1form");
   TF2  *f2form = (TF2*)f.Get("f2form");
   TH1F *h1form = (TH1F*)f.Get("h1form");
   TH2F *h2form = (TH2F*)f.Get("h2form");

   //Divide the canvas in subpads. Plot with different options
   c1->Divide(2,2);
   c1->cd(1);
   f1form->Draw();
   c1->cd(2);
   h1form->Draw();
   c1->cd(3);
   h2form->Draw("box");
   f2form->Draw("cont1same");
   c1->cd(4);
   f2form->Draw("surf");

   ps.Close();

   //count number of lines in ps file
   FILE *fp = fopen("stress.ps","r");
   char line[260];
   Int_t nlines = 0;
   Int_t nlinesGood = 632;
   while (fgets(line,255,fp)) {
      nlines++;
   }
   fclose(fp);
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();
   Bool_t OK = kTRUE;
   if (nlines < nlinesGood-110 || nlines > nlinesGood+110) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s nlines in stress.ps file = %d nlinesGood =%d\n"," ",nlines,nlinesGood);
   }
   delete c1;
   if (gPrintSubBench) { printf("Test  5 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
   
}

//_______________________________________________________________
void stress6()
{
// Test subdirectories in a Root file
// Create many TH1S histograms, make operations between them

   Bprint(6,"Test subdirectories in a Root file");

   TFile f("stress.root","update");
   // create a new subdirectory for each plane
   gRandom->SetSeed(65539);
   const Int_t nplanes = 10;
   const Int_t ncounters = 100;
   char dirname[50];
   char hname[20];
   char htitle[80];
   TH1S *hn[ncounters];
   TH1S *hs[ncounters];
   Int_t i,j,k,id;
   TH1F *hsumPlanes = new TH1F("hsumPlanes","Sum of all planes",100,0,100);
   //Create a subdirectory per detector plane
   for (i=0;i<nplanes;i++) {
      snprintf(dirname,50,"plane%d",i);
      TDirectory *cdplane = f.mkdir(dirname);
      if (cdplane == 0) continue;
      cdplane->cd();
      // create counter histograms
      for (j=0;j<ncounters;j++) {
         snprintf(hname,20,"h%d_%dN",i,j);
         snprintf(htitle,80,"hist for counter:%d in plane:%d North",j,i);
         hn[j] = new TH1S(hname,htitle,100,0,100);
         snprintf(hname,20,"h%d_%dS",i,j);
         snprintf(htitle,80,"hist for counter:%d in plane:%d South",j,i);
         hs[j] = new TH1S(hname,htitle,100,0,100);
      }
      // fill counter histograms randomly
      for (k=0;k<10000;k++) {
         id = Int_t(ncounters*gRandom->Rndm());
         hn[id]->Fill(gRandom->Gaus(60,10));
         hs[id]->Fill(gRandom->Gaus(40,5));
      }
      // Write all objects in directory in memory to disk
      cdplane->Write();
      // Delete all objects from memory
      cdplane->GetList()->Delete();
      f.cd();
   }
   // Now read back all objects from all subdirectories
   // Add North and south histograms in hsumPlanes
   for (i=0;i<nplanes;i++) {
      snprintf(dirname,50,"plane%d",i);
      f.cd(dirname);
      for (j=0;j<ncounters;j++) {
         snprintf(hname,20,"h%d_%dN",i,j);
         TH1S *hnorth; gDirectory->GetObject(hname,hnorth);
         snprintf(hname,20,"h%d_%dS",i,j);
         TH1S *hsouth; gDirectory->GetObject(hname,hsouth);
         if (hnorth == 0 || hsouth == 0) continue;
         hsumPlanes->Add(hnorth);
         hsumPlanes->Add(hsouth);
         delete hnorth; delete hsouth;
      }
      f.cd();    // change current directory to top
   }
   // Verify number of entries, rms and mean value
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();
   Int_t nentries = (Int_t)hsumPlanes->GetEntries();
   Double_t rms   = hsumPlanes->GetRMS();
   Double_t mean  = hsumPlanes->GetMean();
   Int_t nentriesGood = 200000;
   Double_t rmsGood  = 12.745;
   Double_t meanGood = 50.01;
   Double_t diffrms  = TMath::Abs(rmsGood -rms)/rmsGood;
   Double_t diffmean = TMath::Abs(meanGood -mean)/meanGood;
   Bool_t OK = kTRUE;
   if (nentriesGood != nentries || diffrms > 1.e-2 || diffmean > 1.e-2) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s nentries=%d, diffmean=%g, diffrms=%g\n"," ",nentries,diffmean,diffrms);
   }
   if (gPrintSubBench) { printf("Test  6 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress7()
{
// Test TNtuple class with several selection mechanisms
// Test expression cuts
// Test graphical cuts
// Test event lists and operations on event lists
// Compare results of TTree::Draw with results of an explict loop

   Bprint(7,"TNtuple, selections, TCut, TCutG, TEventList");

   TFile f("stress.root","update");
   // Create and fill a TNtuple
   gRandom->SetSeed(65539);
   TNtuple *ntuple = new TNtuple("ntuple","Demo ntuple","px:py:pz:random:i");
   Float_t px, py, pz;
   Int_t nall = 50000;
   Int_t i;
   for (i = 0; i < nall; i++) {
      gRandom->Rannor(px,py);
      pz = px*px + py*py;
      Float_t random = gRandom->Rndm(1);
      ntuple->Fill(px,py,pz,random,i);
   }
   ntuple->Write();

   // Create a graphical cut. Select only events in cut
   TCutG *cutg = new TCutG("cutg",9);
   cutg->SetVarX("py");
   cutg->SetVarY("px");
   cutg->SetPoint(0,-1.75713,2.46193);
   cutg->SetPoint(1,-2.58656,-0.786802);
   cutg->SetPoint(2,-0.179195,-0.101523);
   cutg->SetPoint(3,2.12702,-1.49746);
   cutg->SetPoint(4,2.2484,1.95431);
   cutg->SetPoint(5,0.630004,0.583756);
   cutg->SetPoint(6,-0.381495,2.28426);
   cutg->SetPoint(7,-1.27161,1.01523);
   cutg->SetPoint(8,-1.75713,2.46193);
   TH2F *hpxpy = new TH2F("hpxpy","px vx py with cutg",40,-4,4,40,-4,4);
   ntuple->Draw("px:py>>hpxpy","cutg","goff");
   Int_t npxpy = (Int_t)hpxpy->GetEntries();
   Int_t npxpyGood = 27918;
   hpxpy->Write();
   cutg->Write();
   delete cutg;

   // Fill a TEventList using the standard cut
   ntuple->Draw(">>elist","py<0 && pz>4 && random<0.5","goff");
   TEventList *elist; gDirectory->GetObject("elist",elist);
   // Fill hist htemp using the standard cut
   ntuple->Draw("px>>htemp0","py<0 && pz>4 && random<0.5","goff");
   TH1F *htemp0;  gDirectory->GetObject("htemp0",htemp0);
   Double_t pxmean0 = htemp0->GetMean();
   Double_t pxrms0  = htemp0->GetRMS();

   // Fill hist hcut using a TCut = the standard cut
   TCut cut1 = "py<0 && pz>4 && random<0.5";
   TCut vcut = "px>>hcut";
   ntuple->Draw(vcut,cut1,"goff");
   // Fill hist helist looping on the eventlist in TTree::Draw
   ntuple->SetEventList(elist);
   ntuple->Draw("px>>helist","","goff");
   ntuple->SetEventList(0);
   TH1F *hcut;   gDirectory->GetObject("hcut",hcut);
   TH1F *helist; gDirectory->GetObject("helist",helist);
   Int_t n1 = (Int_t)hcut->GetEntries();
   Int_t n2 = (Int_t)helist->GetEntries();
   htemp0->Write();
   cut1.Write();
   helist->Write();
   hcut->Write();

   // now loop on eventlist explicitly and fill helist again
   Float_t pxr;
   ntuple->SetBranchAddress("px",&pxr);
   TH1F *helistc = (TH1F*)helist->Clone();
   helistc->Reset();
   helistc->SetName("helistc");
   Int_t nlist = elist->GetN();
   for (i=0;i<nlist;i++) {
      Long64_t event = elist->GetEntry(i);
      ntuple->GetEntry(event);
      helistc->Fill(pxr);
   }
   Int_t n3 = (Int_t)helistc->GetEntries();
   Double_t pxmean2 = helistc->GetMean();
   Double_t pxrms2  = helistc->GetRMS();
   helistc->Write();
   elist->Write();

   // Generate several TEventlist objects + total and save them
   char elistname[20];
   char cutname[20];
   TEventList *el[10];
   TEventList *elistall = new TEventList("elistall","Sum of all cuts");
   for (i=0;i<10;i++) {
      snprintf(elistname,20,">>elist%d",i);
      snprintf(cutname,20,"i 10 == %d",i); cutname[1] ='%';
      ntuple->Draw(elistname,cutname,"goff");
      gDirectory->GetObject(&elistname[2],el[i]);
      el[i]->Write();
      elistall->Add(el[i]);
   }
   elistall->Write();

   // Read big list from file and check that the distribution with the list
   // correspond to all events (no cuts)
   delete ntuple;
   TNtuple *nt; gDirectory->GetObject("ntuple",nt);
   nt->SetBranchAddress("px",&pxr);
   TH1F *hpx = new TH1F("hpx","hpx",100,-3,3);
   nt->Draw("px>>hpx","","goff");
   TEventList *all; gDirectory->GetObject("elistall",all);
   nt->SetEstimate(nall); //must be done because the order in eventlist is different
   nt->SetEventList(all);
   TH1F *hall = (TH1F*)hpx->Clone();
   hall->SetName("hall");
   nt->Draw("px>>hall","","goff");
   // Take the difference between the two histograms. Must be empty
   //TH1F hcomp = (*hall) - (*hpx);
   //Double_t compsum = hcomp.GetSum();
   hall->Add(hpx,-1);
   Double_t compsum = hall->GetSum();
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();

   // We can compare entries, means and rms
   Bool_t OK = kTRUE;
   if (n1 != n2 || n1 != n3 || n3 != nlist || nall !=elistall->GetN()
                || npxpy != npxpyGood
                || compsum != 0
                || TMath::Abs(pxmean0-pxmean2) > 0.1
                || TMath::Abs(pxrms0-pxrms2) > 0.01) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s n1=%d, n2=%d, n3=%d, elistallN=%d\n"," ",n1,n2,n3,elistall->GetN());
      printf("%-8s pxmean0=%g, pxmean2=%g, pxrms0=%g\n"," ",pxmean0,pxmean2,pxrms0);
      printf("%-8s pxrms2=%g, compsum=%g, npxpy=%d\n"," ",pxrms2,compsum,npxpy);
   }
   if (gPrintSubBench) { printf("Test  7 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
Int_t stress8read(Int_t nevent)
{
//  Read the event file
//  Loop on all events in the file (reading everything).
//  Count number of bytes read

   TFile *hfile = new TFile("Event.root");
   TTree *tree; hfile->GetObject("T",tree);
   Event *event = 0;
   tree->SetBranchAddress("event",&event);
   Int_t nentries = (Int_t)tree->GetEntries();
   Int_t nev = TMath::Max(nevent,nentries);
   //activate the treeCache
   Int_t cachesize = 10000000; //this is the default value: 10 MBytes
   tree->SetCacheSize(cachesize);
   TTreeCache::SetLearnEntries(1); //one entry is sufficient to learn
   TTreeCache *tc = (TTreeCache*)hfile->GetCacheRead();
   tc->SetEntryRange(0,nevent);
   Int_t nb = 0;
   for (Int_t ev = 0; ev < nev; ev++) {
      nb += tree->GetEntry(ev);        //read complete event in memory
   }
   ntotin  += hfile->GetBytesRead();

   delete event;
   delete hfile;
   return nb;
}


//_______________________________________________________________
Int_t stress8write(Int_t nevent, Int_t comp, Int_t split)
{
//  Create the Event file in various modes
   // comp = compression level
   // split = 1 split mode, 0 = no split

   // Create the Event file, the Tree and the branches
   TFile *hfile = new TFile("Event.root","RECREATE","TTree benchmark ROOT file");
   hfile->SetCompressionLevel(comp);

   // Create one event
   Event *event = new Event();

   // Create a ROOT Tree and one superbranch
   TTree *tree = new TTree("T","An example of a ROOT tree");
   tree->SetAutoSave(100000000);  // autosave when 100 Mbytes written
   Int_t bufsize = 64000;
   if (split)  bufsize /= 4;
   tree->Branch("event", &event, bufsize,split);

   //Fill the Tree
   Int_t ev, nb=0, meanTracks=600;
   Float_t ptmin = 1;
   for (ev = 0; ev < nevent; ev++) {
      event->Build(ev,meanTracks,ptmin);

      nb += tree->Fill();  //fill the tree
   }
   hfile->Write();
   ntotout += hfile->GetBytesWritten();
   delete event;
   delete hfile;
   return nb;
}


//_______________________________________________________________
void stress8(Int_t nevent)
{
//  Run the $ROOTSYS/test/Event program in several configurations.

   Bprint(8,"Trees split and compression modes");

  // First step: make sure the Event shared library exists
  // This test dynamic linking when running in interpreted mode
   if (!TClassTable::GetDict("Event")) {
      Int_t st1 = -1;
      if (gSystem->DynamicPathName("$ROOTSYS/test/libEvent",kTRUE)) {
         st1 = gSystem->Load("$(ROOTSYS)/test/libEvent");
      }
      if (st1 == -1) {
         if (gSystem->DynamicPathName("test/libEvent",kTRUE)) {
            st1 = gSystem->Load("test/libEvent");
         }
         if (st1 == -1) {
            printf("===>stress8 will try to build the libEvent library\n");
            Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
            if (UNIX) gSystem->Exec("(cd $ROOTSYS/test; make Event)");
            else      gSystem->Exec("(cd %ROOTSYS%\\test && nmake libEvent.dll)");
            st1 = gSystem->Load("$(ROOTSYS)/test/libEvent");
         }
      }
   }

   // Create the file not compressed, in no-split mode and read it back
   gRandom->SetSeed(65539);
   Int_t nbw0 = stress8write(100,0,0);
   Int_t nbr0 = stress8read(0);
   Event::Reset();

   // Create the file compressed, in no-split mode and read it back
   gRandom->SetSeed(65539);
   Int_t nbw1 = stress8write(100,1,0);
   Int_t nbr1 = stress8read(0);
   Event::Reset();

   // Create the file compressed, in split mode and read it back
   gRandom->SetSeed(65539);
   Int_t nbw2 = stress8write(nevent,1,9);
   Int_t nbr2 = stress8read(0);
   Event::Reset();

   Bool_t OK = kTRUE;
   if (nbw0 != nbr0 || nbw1 != nbr1 || nbw2 != nbr2) OK = kFALSE;
   if (nbw0 != nbw1) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s nbw0=%d, nbr0=%d, nbw1=%d\n"," ",nbw0,nbr0,nbw1);
      printf("%-8s nbr1=%d, nbw2=%d, nbr2=%d\n"," ",nbr1,nbw2,nbr2);
   }
   if (gPrintSubBench) { printf("Test  8 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
Int_t HistCompare(TH1 *h1, TH1 *h2)
{
// Compare histograms h1 and h2
// Check number of entries, mean and rms
// if means differ by more than 1/1000 of the range return -1
// if rms differs in percent by more than 1/1000 return -2
// Otherwise return difference of number of entries

   Int_t n1       = (Int_t)h1->GetEntries();
   Double_t mean1 = h1->GetMean();
   Double_t rms1  = h1->GetRMS();
   Int_t n2       = (Int_t)h2->GetEntries();
   Double_t mean2 = h2->GetMean();
   Double_t rms2  = h2->GetRMS();
   Float_t xrange = h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin();
   if (TMath::Abs((mean1-mean2)/xrange) > 0.001*xrange) return -1;
   if (rms1 && TMath::Abs((rms1-rms2)/rms1) > 0.001)    return -2;
   return n1-n2;
}

//_______________________________________________________________
void stress9tree(TTree *tree, Int_t realTestNum)
{
// Test selections via TreeFormula
// tree is a TTree when called by stress9
// tree is a TChain when called from stress11
// This is a quite complex test checking the results of TTree::Draw
// or TChain::Draw with an explicit loop on events.
// Also a good test for the interpreter

   Event *event = 0;
   tree->SetBranchAddress("event",&event);
   gROOT->cd();
   TDirectory *hfile = gDirectory;
   Double_t nrsave = TFile::GetFileBytesRead();

   // Each tree->Draw generates an histogram
   tree->Draw("fNtrack>>hNtrack",    "","goff");
   tree->Draw("fNseg>>hNseg",        "","goff");
   tree->Draw("fTemperature>>hTemp", "","goff");
   tree->Draw("fH.GetMean()>>hHmean","","goff");
   tree->Draw("fTracks.fPx>>hPx","fEvtHdr.fEvtNum%10 == 0","goff");
   tree->Draw("fTracks.fPy>>hPy","fEvtHdr.fEvtNum%10 == 0","goff");
   tree->Draw("fTracks.fPz>>hPz","fEvtHdr.fEvtNum%10 == 0","goff");
   tree->Draw("fRandom>>hRandom","fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fMass2>>hMass2",  "fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fBx>>hBx",        "fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fBy>>hBy",        "fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fXfirst>>hXfirst","fEvtHdr.fEvtNum%10 == 2","goff");
   tree->Draw("fYfirst>>hYfirst","fEvtHdr.fEvtNum%10 == 2","goff");
   tree->Draw("fZfirst>>hZfirst","fEvtHdr.fEvtNum%10 == 2","goff");
   tree->Draw("fXlast>>hXlast",  "fEvtHdr.fEvtNum%10 == 3","goff");
   tree->Draw("fYlast>>hYlast",  "fEvtHdr.fEvtNum%10 == 3","goff");
   tree->Draw("fZlast>>hZlast",  "fEvtHdr.fEvtNum%10 == 3","goff");
   tree->Draw("fCharge>>hCharge","fPx < 0","goff");
   tree->Draw("fNpoint>>hNpoint","fPx < 0","goff");
   tree->Draw("fValid>>hValid",  "fPx < 0","goff");

   tree->Draw("fMatrix>>hFullMatrix","","goff");
   tree->Draw("fMatrix[][0]>>hColMatrix","","goff");
   tree->Draw("fMatrix[1][]>>hRowMatrix","","goff");
   tree->Draw("fMatrix[2][2]>>hCellMatrix","","goff");

   tree->Draw("fMatrix - fVertex>>hFullOper","","goff");
   tree->Draw("fMatrix[2][1] - fVertex[5][1]>>hCellOper","","goff");
   tree->Draw("fMatrix[][1]  - fVertex[5][1]>>hColOper","","goff");
   tree->Draw("fMatrix[2][]  - fVertex[5][2]>>hRowOper","","goff");
   tree->Draw("fMatrix[2][]  - fVertex[5][]>>hMatchRowOper","","goff");
   tree->Draw("fMatrix[][2]  - fVertex[][1]>>hMatchColOper","","goff");
   tree->Draw("fMatrix[][2]  - fVertex[][]>>hRowMatOper","","goff");
   tree->Draw("fMatrix[][2]  - fVertex[5][]>>hMatchDiffOper","","goff");
   tree->Draw("fMatrix[][]   - fVertex[][]>>hFullOper2","","goff");

   if (gPrintSubBench) { printf("\n"); printf("Test %2dD: ",realTestNum); gBenchmark->Show("stress");gBenchmark->Start("stress"); }

   ntotin  += TFile::GetFileBytesRead() -nrsave;

   //Get pointers to the histograms generated above
   TH1F *hNtrack = (TH1F*)hfile->Get("hNtrack");
   TH1F *hNseg   = (TH1F*)hfile->Get("hNseg");
   TH1F *hTemp   = (TH1F*)hfile->Get("hTemp");
   TH1F *hHmean  = (TH1F*)hfile->Get("hHmean");
   TH1F *hPx     = (TH1F*)hfile->Get("hPx");
   TH1F *hPy     = (TH1F*)hfile->Get("hPy");
   TH1F *hPz     = (TH1F*)hfile->Get("hPz");
   TH1F *hRandom = (TH1F*)hfile->Get("hRandom");
   TH1F *hMass2  = (TH1F*)hfile->Get("hMass2");
   TH1F *hBx     = (TH1F*)hfile->Get("hBx");
   TH1F *hBy     = (TH1F*)hfile->Get("hBy");
   TH1F *hXfirst = (TH1F*)hfile->Get("hXfirst");
   TH1F *hYfirst = (TH1F*)hfile->Get("hYfirst");
   TH1F *hZfirst = (TH1F*)hfile->Get("hZfirst");
   TH1F *hXlast  = (TH1F*)hfile->Get("hXlast");
   TH1F *hYlast  = (TH1F*)hfile->Get("hYlast");
   TH1F *hZlast  = (TH1F*)hfile->Get("hZlast");
   TH1F *hCharge = (TH1F*)hfile->Get("hCharge");
   TH1F *hNpoint = (TH1F*)hfile->Get("hNpoint");
   TH1F *hValid  = (TH1F*)hfile->Get("hValid");

   TH1F *hFullMatrix    = (TH1F*)hfile->Get("hFullMatrix");
   TH1F *hColMatrix     = (TH1F*)hfile->Get("hColMatrix");
   TH1F *hRowMatrix     = (TH1F*)hfile->Get("hRowMatrix");
   TH1F *hCellMatrix    = (TH1F*)hfile->Get("hCellMatrix");
   TH1F *hFullOper      = (TH1F*)hfile->Get("hFullOper");
   TH1F *hCellOper      = (TH1F*)hfile->Get("hCellOper");
   TH1F *hColOper       = (TH1F*)hfile->Get("hColOper");
   TH1F *hRowOper       = (TH1F*)hfile->Get("hRowOper");
   TH1F *hMatchRowOper  = (TH1F*)hfile->Get("hMatchRowOper");
   TH1F *hMatchColOper  = (TH1F*)hfile->Get("hMatchColOper");
   TH1F *hRowMatOper    = (TH1F*)hfile->Get("hRowMatOper");
   TH1F *hMatchDiffOper = (TH1F*)hfile->Get("hMatchDiffOper");
   TH1F *hFullOper2     = (TH1F*)hfile->Get("hFullOper2");

   //We make clones of the generated histograms
   //We set new names and reset the clones.
   //We want to have identical histogram limits
   TH1F *bNtrack = (TH1F*)hNtrack->Clone(); bNtrack->SetName("bNtrack"); bNtrack->Reset();
   TH1F *bNseg   = (TH1F*)hNseg->Clone();   bNseg->SetName("bNseg");     bNseg->Reset();
   TH1F *bTemp   = (TH1F*)hTemp->Clone();   bTemp->SetName("bTemp");     bTemp->Reset();
   TH1F *bHmean  = (TH1F*)hHmean->Clone();  bHmean->SetName("bHmean");   bHmean->Reset();
   TH1F *bPx     = (TH1F*)hPx->Clone();     bPx->SetName("bPx");         bPx->Reset();
   TH1F *bPy     = (TH1F*)hPy->Clone();     bPy->SetName("bPy");         bPy->Reset();
   TH1F *bPz     = (TH1F*)hPz->Clone();     bPz->SetName("bPz");         bPz->Reset();
   TH1F *bRandom = (TH1F*)hRandom->Clone(); bRandom->SetName("bRandom"); bRandom->Reset();
   TH1F *bMass2  = (TH1F*)hMass2->Clone();  bMass2->SetName("bMass2");   bMass2->Reset();
   TH1F *bBx     = (TH1F*)hBx->Clone();     bBx->SetName("bBx");         bBx->Reset();
   TH1F *bBy     = (TH1F*)hBy->Clone();     bBy->SetName("bBy");         bBy->Reset();
   TH1F *bXfirst = (TH1F*)hXfirst->Clone(); bXfirst->SetName("bXfirst"); bXfirst->Reset();
   TH1F *bYfirst = (TH1F*)hYfirst->Clone(); bYfirst->SetName("bYfirst"); bYfirst->Reset();
   TH1F *bZfirst = (TH1F*)hZfirst->Clone(); bZfirst->SetName("bZfirst"); bZfirst->Reset();
   TH1F *bXlast  = (TH1F*)hXlast->Clone();  bXlast->SetName("bXlast");   bXlast->Reset();
   TH1F *bYlast  = (TH1F*)hYlast->Clone();  bYlast->SetName("bYlast");   bYlast->Reset();
   TH1F *bZlast  = (TH1F*)hZlast->Clone();  bZlast->SetName("bZlast");   bZlast->Reset();
   TH1F *bCharge = (TH1F*)hCharge->Clone(); bCharge->SetName("bCharge"); bCharge->Reset();
   TH1F *bNpoint = (TH1F*)hNpoint->Clone(); bNpoint->SetName("bNpoint"); bNpoint->Reset();
   TH1F *bValid  = (TH1F*)hValid->Clone();  bValid->SetName("bValid");   bValid->Reset();

   TH1F *bFullMatrix    =(TH1F*)hFullMatrix->Clone();    bFullMatrix->SetName("bFullMatrix");       bFullMatrix->Reset();
   TH1F *bColMatrix    = (TH1F*)hColMatrix->Clone();     bColMatrix->SetName("bColMatrix");         bColMatrix->Reset();
   TH1F *bRowMatrix    = (TH1F*)hRowMatrix->Clone();     bRowMatrix->SetName("bRowMatrix");         bRowMatrix->Reset();
   TH1F *bCellMatrix   = (TH1F*)hCellMatrix->Clone();    bCellMatrix->SetName("bCellMatrix");       bCellMatrix->Reset();
   TH1F *bFullOper     = (TH1F*)hFullOper->Clone();      bFullOper->SetName("bFullOper");           bFullOper->Reset();
   TH1F *bCellOper     = (TH1F*)hCellOper->Clone();      bCellOper->SetName("bCellOper");           bCellOper->Reset();
   TH1F *bColOper      = (TH1F*)hColOper->Clone();       bColOper->SetName("bColOper");             bColOper->Reset();
   TH1F *bRowOper      = (TH1F*)hRowOper->Clone();       bRowOper->SetName("bRowOper");             bRowOper->Reset();
   TH1F *bMatchRowOper = (TH1F*)hMatchRowOper->Clone();  bMatchRowOper->SetName("bMatchRowOper");   bMatchRowOper->Reset();
   TH1F *bMatchColOper = (TH1F*)hMatchColOper->Clone();  bMatchColOper->SetName("bMatchColOper");   bMatchColOper->Reset();
   TH1F *bRowMatOper   = (TH1F*)hRowMatOper->Clone();    bRowMatOper->SetName("bRowMatOper");       bRowMatOper->Reset();
   TH1F *bMatchDiffOper= (TH1F*)hMatchDiffOper->Clone(); bMatchDiffOper->SetName("bMatchDiffOper"); bMatchDiffOper->Reset();
   TH1F *bFullOper2    = (TH1F*)hFullOper2->Clone();     bFullOper2->SetName("bFullOper2");         bFullOper2->Reset();

   // Loop with user code on all events and fill the b histograms
   // The code below should produce identical results to the tree->Draw above

   TClonesArray *tracks = event->GetTracks();
   Int_t nev = (Int_t)tree->GetEntries();
   Int_t i, ntracks, evmod,i0,i1;
   Track *t;
   EventHeader *head;
   Int_t nbin = 0;
   for (Int_t ev=0;ev<nev;ev++) {
      nbin += tree->GetEntry(ev);
      head = event->GetHeader();
      evmod = head->GetEvtNum()%10;
      bNtrack->Fill(event->GetNtrack());
      bNseg->Fill(event->GetNseg());
      bTemp->Fill(event->GetTemperature());
      bHmean->Fill(event->GetHistogram()->GetMean());
      ntracks = event->GetNtrack();
      for(i0=0;i0<4;i0++) {
         for(i1=0;i1<4;i1++) {
            bFullMatrix->Fill(event->GetMatrix(i0,i1));
         }
         bColMatrix->Fill(event->GetMatrix(i0,0));
         bRowMatrix->Fill(event->GetMatrix(1,i0)); // done here because the matrix is square!
      }
      bCellMatrix->Fill(event->GetMatrix(2,2));
      if ( 5 < ntracks ) {
         t = (Track*)tracks->UncheckedAt(5);
         for(i0=0;i0<4;i0++) {
            for(i1=0;i1<4;i1++) {
            }
            bColOper->Fill( event->GetMatrix(i0,1) - t->GetVertex(1) );
            bRowOper->Fill( event->GetMatrix(2,i0) - t->GetVertex(2) );
         }
         for(i0=0;i0<3;i0++) {
            bMatchRowOper->Fill( event->GetMatrix(2,i0) - t->GetVertex(i0) );
            bMatchDiffOper->Fill( event->GetMatrix(i0,2) - t->GetVertex(i0) );
         }
         bCellOper->Fill( event->GetMatrix(2,1) - t->GetVertex(1) );
      }
      for (i=0;i<ntracks;i++) {
         t = (Track*)tracks->UncheckedAt(i);
         if (evmod == 0) bPx->Fill(t->GetPx());
         if (evmod == 0) bPy->Fill(t->GetPy());
         if (evmod == 0) bPz->Fill(t->GetPz());
         if (evmod == 1) bRandom->Fill(t->GetRandom());
         if (evmod == 1) bMass2->Fill(t->GetMass2());
         if (evmod == 1) bBx->Fill(t->GetBx());
         if (evmod == 1) bBy->Fill(t->GetBy());
         if (evmod == 2) bXfirst->Fill(t->GetXfirst());
         if (evmod == 2) bYfirst->Fill(t->GetYfirst());
         if (evmod == 2) bZfirst->Fill(t->GetZfirst());
         if (evmod == 3) bXlast->Fill(t->GetXlast());
         if (evmod == 3) bYlast->Fill(t->GetYlast());
         if (evmod == 3) bZlast->Fill(t->GetZlast());
         if (t->GetPx() < 0) {
            bCharge->Fill(t->GetCharge());
            bNpoint->Fill(t->GetNpoint());
            bValid->Fill(t->GetValid());
         }
         if (i<4) {
            for(i1=0;i1<3;i1++) { // 3 is the min of the 2nd dim of Matrix and Vertex
               bFullOper ->Fill( event->GetMatrix(i,i1) - t->GetVertex(i1) );
               bFullOper2->Fill( event->GetMatrix(i,i1) - t->GetVertex(i1) );
               bRowMatOper->Fill( event->GetMatrix(i,2) - t->GetVertex(i1) );
            }
            bMatchColOper->Fill( event->GetMatrix(i,2) - t->GetVertex(1) );
         }
      }
   }

   // Compare h and b histograms
   Int_t cNtrack = HistCompare(hNtrack,bNtrack);
   Int_t cNseg   = HistCompare(hNseg,bNseg);
   Int_t cTemp   = HistCompare(hTemp,bTemp);
   Int_t cHmean  = HistCompare(hHmean,bHmean);
   Int_t cPx     = HistCompare(hPx,bPx);
   Int_t cPy     = HistCompare(hPy,bPy);
   Int_t cPz     = HistCompare(hPz,bPz);
   Int_t cRandom = HistCompare(hRandom,bRandom);
   Int_t cMass2  = HistCompare(hMass2,bMass2);
   Int_t cBx     = HistCompare(hBx,bBx);
   Int_t cBy     = HistCompare(hBy,bBy);
   Int_t cXfirst = HistCompare(hXfirst,bXfirst);
   Int_t cYfirst = HistCompare(hYfirst,bYfirst);
   Int_t cZfirst = HistCompare(hZfirst,bZfirst);
   Int_t cXlast  = HistCompare(hXlast,bXlast);
   Int_t cYlast  = HistCompare(hYlast,bYlast);
   Int_t cZlast  = HistCompare(hZlast,bZlast);
   Int_t cCharge = HistCompare(hCharge,bCharge);
   Int_t cNpoint = HistCompare(hNpoint,bNpoint);
   Int_t cValid  = HistCompare(hValid,bValid);

   Int_t cFullMatrix   = HistCompare(hFullMatrix,bFullMatrix);
   Int_t cColMatrix    = HistCompare(hColMatrix,bColMatrix);
   Int_t cRowMatrix    = HistCompare(hRowMatrix,bRowMatrix);
   Int_t cCellMatrix   = HistCompare(hCellMatrix,bCellMatrix);
   Int_t cFullOper     = HistCompare(hFullOper,bFullOper);
   Int_t cCellOper     = HistCompare(hCellOper,bCellOper);
   Int_t cColOper      = HistCompare(hColOper,bColOper);
   Int_t cRowOper      = HistCompare(hRowOper,bRowOper);
   Int_t cMatchRowOper = HistCompare(hMatchRowOper,bMatchRowOper);
   Int_t cMatchColOper = HistCompare(hMatchColOper,bMatchColOper);
   Int_t cRowMatOper   = HistCompare(hRowMatOper,bRowMatOper);
   Int_t cMatchDiffOper= HistCompare(hMatchDiffOper,bMatchDiffOper);
   Int_t cFullOper2    = HistCompare(hFullOper2,bFullOper2);

   delete event;
   Event::Reset();
   ntotin += nbin;

   if (gPrintSubBench) { 
      printf("Test %2dC: ",realTestNum); 
      gBenchmark->Show("stress");gBenchmark->Start("stress");
      // Since we disturbed the flow (due to the double benchmark printing),
      // let's repeat the header!
      printf("Test %2d : ",realTestNum);
   }
   
   Bool_t OK = kTRUE;
   if (cNtrack || cNseg   || cTemp  || cHmean || cPx    || cPy     || cPz) OK = kFALSE;
   if (cRandom || cMass2  || cBx    || cBy    || cXfirst|| cYfirst || cZfirst) OK = kFALSE;
   if (cXlast  || cYlast  || cZlast || cCharge|| cNpoint|| cValid) OK = kFALSE;
   if (cFullMatrix || cColMatrix || cRowMatrix || cCellMatrix || cFullOper ) OK = kFALSE;
   if (cCellOper || cColOper || cRowOper || cMatchRowOper || cMatchColOper ) OK = kFALSE;
   if (cRowMatOper || cMatchDiffOper || cFullOper2 ) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s cNtrak =%d, cNseg  =%d, cTemp  =%d, cHmean =%d\n"," ",cNtrack,cNseg,cTemp,cHmean);
      printf("%-8s cPx    =%d, cPy    =%d, cPz    =%d, cRandom=%d\n"," ",cPx,cPy,cPz,cRandom);
      printf("%-8s cMass2 =%d, cbx    =%d, cBy    =%d, cXfirst=%d\n"," ",cMass2,cBx,cBy,cXfirst);
      printf("%-8s cYfirst=%d, cZfirst=%d, cXlast =%d, cYlast =%d\n"," ",cYfirst,cZfirst,cXlast,cYlast);
      printf("%-8s cZlast =%d, cCharge=%d, cNpoint=%d, cValid =%d\n"," ",cZlast,cCharge,cNpoint,cValid);
      printf("%-8s cFullMatrix=%d, cColMatrix=%d, cRowMatrix=%d, cCellMatrix=%d\n"," ",cFullMatrix,cColMatrix,cRowMatrix,cCellMatrix);
      printf("%-8s cFullOper=%d, cCellOper=%d, cColOper=%d, cRowOper=%d\n"," ",cFullOper,cCellOper,cColOper,cRowOper);
      printf("%-8s cMatchRowOper=%d, cMatchColOper=%d, cRowMatOper=%d, cMatchDiffOper=%d\n"," ",cMatchRowOper,cMatchColOper,cRowMatOper,cMatchDiffOper);
      printf("%-8s cFullOper2=%d\n"," ",cFullOper2);
   }
}

//_______________________________________________________________
void stress9()
{
// Analyse the file Event.root generated in the last part of test8

   Bprint(9,"Analyze Event.root file of stress 8");

   gROOT->GetList()->Delete();
   TFile *hfile = new TFile("Event.root");
   TTree *tree; hfile->GetObject("T",tree);

   stress9tree(tree,9);

   // Save test9 histograms
   TFile f("stress_test9.root","recreate");
   gROOT->GetList()->Write();
   gROOT->GetList()->Delete();
   ntotout += f.GetBytesWritten();


   delete hfile;
}

//_______________________________________________________________
void stress10()
{
// Make 10 Trees starting from the Event.root tree.
// Events for which event_number%10 == 0 go to Event_0.root
// Events for which event_number%10 == 1 go to Event_1.root
//...
// Events for which event_number%10 == 9 go to Event_9.root

   Bprint(10,"Create 10 files starting from Event.root");

   TFile *hfile = new TFile("Event.root");
   if (hfile==0 || hfile->IsZombie()) {
      delete hfile;
      printf("failed\n");
      return;
   }
   TTree *tree; hfile->GetObject("T",tree);

   Event *event = 0;
   tree->SetBranchAddress("event",&event);

   // Create 10 clones of this tree
   char filename[20];
   TTree *chTree[10];
   TFile *chfile[10];
   Int_t file;
   for (file=0;file<10;file++) {
      snprintf(filename,20,"Event_%d.root",file);
      chfile[file] = new TFile(filename,"recreate");
      if (file>=5) {
         chfile[file]->SetCompressionAlgorithm(ROOT::kLZMA);
      }
      chTree[file] = (TTree*)tree->CloneTree(0);
   }

   // Fill the small trees
   Int_t nev = (Int_t)tree->GetEntries();
   Int_t evmod, nbin=0, nbout=0;
   EventHeader *head;
   for (Int_t ev=0;ev<nev;ev++) {
      nbin += tree->GetEntry(ev);
      head = event->GetHeader();
      evmod = head->GetEvtNum()%10;
      nbout += chTree[evmod]->Fill();
      event->Clear();
   }
   // save headers
   Int_t ntot = 0;
   for (file=0;file<10;file++) {
      ntot += (Int_t)chTree[file]->GetEntries();
      chfile[file]->Write();
      delete chfile[file];
   }
   delete event;
   delete hfile;
   Event::Reset();
   ntotin  += nbin;
   ntotout += nbout;

   //We compare the number of bytes read from the big file
   //with the total number of bytes written in the 10 small files
   Bool_t OK = kTRUE;
   if (nbin != nbout || nev != ntot) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s nbin=%d, nbout=%d, nev=%d, ntot=%d\n"," ",nbin,nbout,nev,ntot);
   }
   if (gPrintSubBench) { printf("Test 10 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress11()
{
// Test chains of Trees
// We make a TChain using the 10 files generated in test10
// We expect the same results when analyzing the chain than
// in the analysis of the original big file Event.root in test9.
// Because TChain derives from TTree, we can use the same
// analysis procedure "stress9tree"

   Bprint(11,"Test chains of Trees using the 10 files");

   gROOT->GetList()->Delete();
   TChain *chain = new TChain("T");
   char filename[20];
   Int_t file;
   for (file=0;file<10;file++) {
      snprintf(filename,20,"Event_%d.root",file);
      chain->Add(filename);
   }

   stress9tree(chain,11);

   // Save test11 histograms
   delete chain;
   TFile f("stress_test11.root","recreate");
   gROOT->GetList()->Write();
   gROOT->GetList()->Delete();
   ntotout += f.GetBytesWritten();
}

//_______________________________________________________________
void stress12(Int_t testid)
{
// Compare histograms of stress9 with stress11

   if (testid == 12) Bprint(12,"Compare histograms of test 9 and 11");

   TFile f9("stress_test9.root");
   TFile f11("stress_test11.root");
   //Let's loop on all keys of second file
   //We expect to find the same keys in the original stress9 file
   TIter next(f11.GetListOfKeys());
   TKey *key;
   TH1F *h9, *h11;
   Int_t comp, ngood = 0;
   while ((key=(TKey*)next())) {
      if (strcmp(key->GetClassName(),"TH1F")) continue; //may be a TList of TStreamerInfo
      h9  = (TH1F*)f9.Get(key->GetName());
      h11 = (TH1F*)f11.Get(key->GetName());
      if (h9 == 0 || h11 == 0) continue;
      comp = HistCompare(h9,h11);
      if (comp == 0) ngood++;
   }
   ntotin += f9.GetBytesRead();
   ntotin += f11.GetBytesRead();
   Bool_t OK = kTRUE;
   if (ngood < 40) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
      printf("%-8s ngood=%d\n"," ",ngood);
   }
   if (gPrintSubBench) { printf("Test 12 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress13()
{
// test of TChain::Merge
// The 10 small Tree files generated in stress10 are again merged
// into one single file.
// Should be the same as the file generated in stress8, except
// that events will be in a different order.
// But global analysis histograms should be identical (checked by stress14)

   Bprint(13,"Test merging files of a chain");

   gROOT->GetList()->Delete();
   TChain *chain = new TChain("T");
   char filename[20];
   Int_t file;
   for (file=0;file<10;file++) {
      snprintf(filename,20,"Event_%d.root",file);
      chain->Add(filename);
   }

   chain->Merge("Event.root");

   Double_t chentries = chain->GetEntries();
   delete chain;

   Event::Reset();
   gROOT->GetList()->Delete();

   TFile f("Event.root");
   TTree *tree = (TTree*)f.Get("T");
   ntotin  += (Double_t)f.GetEND();
   ntotout += (Double_t)f.GetEND();

   Bool_t OK = kTRUE;
   if (chentries != tree->GetEntries()) OK = kFALSE;
   if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
   else    {
      printf("failed\n");
   }
   if (gPrintSubBench) { printf("Test 13 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

//_______________________________________________________________
void stress14()
{
// Verify that stress13 has correctly rebuild the original Event.root

   Bprint(14,"Check correct rebuilt of Event.root in test 13");

   stress12(14);
}

//_______________________________________________________________
void stress15()
{
// Divert some branches to separate files

   Bprint(15,"Divert Tree branches to separate files");

   //Get old file, old tree and set top branch address
   //We want to copy only a few branches.
   TFile *oldfile = new TFile("Event.root");
   if (oldfile->IsZombie()) {
      printf("failed\n");
      return;
   }   
   TTree *oldtree; oldfile->GetObject("T",oldtree);
   Event *event   = 0;
   oldtree->SetBranchAddress("event",&event);
   oldtree->SetBranchStatus("*",0);
   oldtree->SetBranchStatus("event",1);
   oldtree->SetBranchStatus("fNtrack",1);
   oldtree->SetBranchStatus("fNseg",1);
   oldtree->SetBranchStatus("fH",1);


   //Create a new file + a clone of old tree header. Do not copy events
   TFile *newfile = new TFile("stress_small.root","recreate");
   TTree *newtree = oldtree->CloneTree(0);

   //Divert branch fH to a separate file and copy all events
   newtree->GetBranch("fH")->SetFile("stress_fH.root");
   newtree->CopyEntries(oldtree);

   newfile->Write();
   ntotin  += oldfile->GetBytesRead();
   ntotout += newfile->GetBytesWritten();
   delete event;
   delete newfile;
   delete oldfile;
   Event::Reset();
   gROOT->GetList()->Delete();

   // Open small file, histogram fNtrack and fH
   newfile = new TFile("stress_small.root");
   newfile->GetObject("T", newtree);
   newtree->Draw("fNtrack>>hNtrack","","goff");
   newtree->Draw("fH.GetMean()>>hHmean","","goff");
   TH1F *hNtrack; newfile->GetObject("hNtrack",hNtrack);
   TH1F *hHmean; newfile->GetObject("hHmean",hHmean);
   ntotin  += newfile->GetBytesRead();

   // Open old reference file of stress9
   oldfile = new TFile("stress_test9.root");
   if (oldfile->IsZombie()) {
      printf("failed\n");
      return;
   }
   TH1F *bNtrack; oldfile->GetObject("bNtrack",bNtrack);
   TH1F *bHmean;  oldfile->GetObject("bHmean",bHmean);
   Int_t cNtrack = HistCompare(hNtrack,bNtrack);
   Int_t cHmean  = HistCompare(hHmean, bHmean);
   delete newfile;
   delete oldfile;
   Event::Reset();
   gROOT->GetList()->Delete();

   Bool_t OK = kTRUE;
   if (cNtrack || cHmean) OK = kFALSE;
    if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
  else    {
      printf("failed\n");
      printf("%-8s cNtrack=%d, cHmean=%d\n"," ",cNtrack,cHmean);
   }
   if (gPrintSubBench) { printf("Test 15 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

void stress16()
{
// Prototype trigger simulation for the LHCb experiment
// This test nested loops with the interpreter.
// Expected to run fast with the compiler, slow with the interpreter.
// This code is extracted from an original macro by Hans Dijkstra (LHCb)
// The program generates histograms and profile histograms.
// A canvas with subpads containing the results is sent to Postscript.
// We check graphics results by counting the number of lines in the ps file.

   Bprint(16,"CINT test (3 nested loops) with LHCb trigger");

   const int nbuf    = 153;    // buffer size
   const int nlev    = 4;      // number of trigger levels
   const int nstep   = 50000;  // number of steps
   const int itt[4]  = { 1000, 4000, 40000, 400000 }; // time needed per trigger
   const float a[4]  = { 0.25, 0.04, 0.25, 0 };       // acceptance/trigger (last always 0)

   int   i, il, istep, itim[192], itrig[192], it, im, ipass;
   float dead, sum[10];

   // create histogram and array of profile histograms
   TCanvas *c = new TCanvas("laten","latency simulation",700,600);
   gROOT->LoadClass("TPostScript","Postscript");
   TPostScript ps("stress_lhcb.ps",112);
   gRandom->SetSeed(65539);
   TFile f("stress_lhcb.root", "recreate");
   TH1F *pipe = new TH1F("pipe", "free in pipeline", nbuf+1, -0.5, nbuf+0.5);
   pipe->SetLineColor(2);
   pipe->SetFillColor(2);

   TProfile *hp[nlev+1];
   TProfile::Approximate();
   for (i = 0; i <= nlev; i++) {
      char s[64];
      snprintf(s,64, "buf%d", i);
      hp[i] = new TProfile(s, "in buffers", 1000, 0,nstep, -1., 1000.);
      hp[i]->SetLineColor(2);
   }

   dead   = 0;
   sum[0] = nbuf;
   for (i = 1; i <= nlev; i++) sum[i] = 0;
   for (i = 0; i < nbuf; i++) { itrig[i] = 0; itim[i] = 0; }

   for (istep = 0; istep < nstep; istep++) {
      // evaluate status of buffer
      pipe->Fill(sum[0]);
      if ((istep+1)%10 == 0) {
         for (i = 0; i <= nlev; i++)
            hp[i]->Fill((float)istep, sum[i], 1.);
      }

      ipass = 0;
      for (i = 0; i < nbuf; i++) {
         it = itrig[i];
         if (it >= 1) {
            // add 25 ns to all times
            itim[i] += 25;
            im = itim[i];
            // level decisions
            for (il = 0; il < nlev; il++) {
               if (it == il+1 && im > itt[il]) {
                  if (gRandom->Rndm() > a[il]) {
                     itrig[i] = -1;
                     sum[0]++;
                     sum[il+1]--;
                  } else {
                     itrig[i]++;
                     sum[il+1]--;
                     sum[il+2]++;
                  }
               }
            }
         } else if (ipass == 0) {
            itrig[i] = 1;
            itim[i]  = 25;
            sum[0]--;
            sum[1]++;
            ipass++;
         }
      }
      if (ipass == 0) dead++;
   }
//   Float_t deadTime = 100.*dead/nstep;

   // View results in the canvas and make the Postscript file

   c->Divide(2,3);
   c->cd(1); pipe->Draw();
   c->cd(2); hp[0]->Draw();
   c->cd(3); hp[1]->Draw();
   c->cd(4); hp[2]->Draw();
   c->cd(5); hp[3]->Draw();
   c->cd(6); hp[4]->Draw();
   ps.Close();

   f.Write();
   ntotout += f.GetBytesWritten();

   // Check length of Postscript file
   FILE *fp = fopen("stress_lhcb.ps","r");
   char line[260];
   Int_t nlines = 0;
   Int_t nlinesGood = 2121;
   Bool_t counting = kFALSE;
   while (fgets(line,255,fp)) {
      if (counting) nlines++;
      if (strstr(line,"%%EndProlog")) counting = kTRUE;
   }
   fclose(fp);
   delete c;
   Bool_t OK = kTRUE;
   if (nlines < nlinesGood-100 || nlines > nlinesGood+100) OK = kFALSE;
    if (OK) 
    {
  #if VERBOSE
      printf("OK\n");
  #endif
    }
  else    {
      printf("failed\n");
      printf("%-8s nlines in stress_lhcb.ps file = %d\n"," ",nlines);
   }
   if (gPrintSubBench) { printf("Test 16 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

void cleanup()
{
  gSystem->Unlink("stress-vmatrix.root");
  gSystem->Unlink("stress-vvector.root");
  gSystem->Unlink("stress-vdecomp.root");

   gSystem->Unlink("Event.root");
   gSystem->Unlink("Event_0.root");
   gSystem->Unlink("Event_1.root");
   gSystem->Unlink("Event_2.root");
   gSystem->Unlink("Event_3.root");
   gSystem->Unlink("Event_4.root");
   gSystem->Unlink("Event_5.root");
   gSystem->Unlink("Event_6.root");
   gSystem->Unlink("Event_7.root");
   gSystem->Unlink("Event_8.root");
   gSystem->Unlink("Event_9.root");
   gSystem->Unlink("stress.ps");
   gSystem->Unlink("stress.root");
   gSystem->Unlink("stress_fH.root");
   gSystem->Unlink("stress_lhcb.ps");
   gSystem->Unlink("stress_lhcb.root");
   gSystem->Unlink("stress_small.root");
   gSystem->Unlink("stress_test9.root");
   gSystem->Unlink("stress_test11.root");

  for  (Int_t i = 0; i < NG; ++i )
  {
    gSystem->Unlink(Form("%s_diff.root",exps[i]));
  }
}

int main(int argc,const char *argv[]) 
{
  // 1st arg = run/round number 
  // 2nd arg = boolean : with or without the i/o stress
  // 3rd arg = path to the directory with the input files

  if (argc<4)
  {
    std::cerr << "Need 3 arguments !" << std::endl;
    std::cerr << "1st arg = run/round number " << std::endl;
    std::cerr << "2nd arg = boolean : with or without the i/o stress" << std::endl;
    std::cerr << "3rd arg = path to the directory with the input files" << std::endl;
    return -1;
  }

  gBenchmark = new TBenchmark();

  std::map<std::string,Double_t> results;

  int n = atoi(argv[1]);
  int withio = atoi(argv[2]);
  
  gInputFiles = argv[3];

  results["Fit"] = stressFit();
  results["Linear"] = stressLinear();  
  results["Spectrum"] = stressSpectrum();
  if ( withio )
  {
    results["Geometry"] = stressGeometry();
    results["Mix"] = stressMix();
  }


  cleanup();

  delete gBenchmark;

  std::map<std::string,Double_t>::const_iterator it;
  Double_t ct(0.0);

  for ( it = results.begin(); it != results.end(); ++it ) 
  {
    ct += it->second;
  }

//  Double_t reftime = 348.3; //pcbrun4 compiled and 490.5 seconds real time
  const Double_t rootmarks = ct; //800*reftime/ct;

  std::cout << Form("%d %7.2f",n,rootmarks) << std::endl;

  return 0;
}

