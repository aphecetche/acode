#ifndef ROOTSTRESS_H
#define ROOTSTRESS_H

#include "Rtypes.h"
#include "TMatrixDfwd.h"

class TTree;

Double_t stressFit(const char *type = "Minuit", const char *algo = "Migrad", Int_t N = 2000);
Double_t stressGeometry(const char *exp="*", Bool_t generate_ref=kFALSE);
Double_t stressLinear(Int_t maxSizeReq=100,Int_t verbose=0);
Double_t stressMix(Int_t nevent=1000, Int_t style=1, Int_t printSubBenchmark=kFALSE, UInt_t portion =65535);
Double_t stressSpectrum(Int_t ntimes=1000);

void StatusPrint(Int_t id,const TString &title, Int_t nsuccess, Int_t nattempts);
void StatusPrint(Int_t id,const TString &title,Bool_t status);

void mstress_allocation            (Int_t msize);
void mstress_matrix_fill           (Int_t rsize,Int_t csize);
void mstress_element_op            (Int_t rsize,Int_t csize);
void mstress_binary_ebe_op         (Int_t rsize, Int_t csize);
void mstress_transposition         (Int_t msize);
void mstress_special_creation      (Int_t dim);
void mstress_matrix_fill           (Int_t rsize,Int_t csize);
void mstress_matrix_promises       (Int_t dim);
void mstress_norms                 (Int_t rsize,Int_t csize);
void mstress_determinant           (Int_t msize);
void mstress_mm_multiplications    ();
void mstress_sym_mm_multiplications(Int_t msize);
void mstress_vm_multiplications    ();
void mstress_inversion             ();
void mstress_matrix_io             ();

void spstress_allocation           (Int_t msize);
void spstress_matrix_fill          (Int_t rsize,Int_t csize);
void spstress_element_op           (Int_t rsize,Int_t csize);
void spstress_binary_ebe_op        (Int_t rsize, Int_t csize);
void spstress_transposition        (Int_t msize);
void spstress_norms                (Int_t rsize,Int_t csize);
void spstress_mm_multiplications   ();
void spstress_vm_multiplications   ();
void spstress_matrix_slices        (Int_t vsize);
void spstress_matrix_io            ();

void vstress_allocation            (Int_t msize);
void vstress_element_op            (Int_t vsize);
void vstress_binary_op             (Int_t vsize);
void vstress_norms                 (Int_t vsize);
void vstress_matrix_slices         (Int_t vsize);
void vstress_vector_io             ();

Bool_t test_svd_expansion          (const TMatrixD &A);
void   astress_decomp              ();
void   astress_lineqn              ();
void   astress_pseudo              ();
void   astress_eigen               (Int_t msize);
void   astress_decomp_io           (Int_t msize);

void   stress_backward_io          ();

void   cleanup                     ();


void stress1();
void stress2();
void stress3();
void stress4();
void stress5();
void stress6();
void stress7();
void stress8(Int_t nevent);
void stress9tree(TTree *tree, Int_t realTestNum);
void stress9();
void stress10();
void stress11();
void stress12(Int_t testid);
void stress13();
void stress14();
void stress15();
void stress16();

void findPeaks(Int_t pmin, Int_t pmax, Int_t &nfound, Int_t &ngood, Int_t &nghost);

#endif
