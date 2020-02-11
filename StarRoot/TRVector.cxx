#include <assert.h>
#include "StarRoot/TRVector.h"
#include "TError.h"
#include "TString.h"


TRVector::TRVector(Int_t nrows):    TRMatrix(nrows, 1) {;}


TRVector::TRVector(Int_t nrows, const Double_t* fArray): TRMatrix(nrows, 1, fArray) {;}


TRVector::TRVector(Int_t nrows, const Float_t* fArray): TRMatrix(nrows, 1, fArray) {;}


TRVector::TRVector(const TRMatrix &A, ETRMatrixCreatorsOp kop, const TRVector &B) :
  TRMatrix(A, kop, B) {}


TRVector::TRVector(const TRMatrix &S, Int_t I)  // get I'th row
{
  if (I == 0) {
    ::Error("TRVector::TRVector(const TRMatrix&)", "index i %d out of bounds (size: %d, this: %p)",
            I, S.NI(), this); I = 1;
  }

  if (I > S.NI()) {
    ::Error("TRVector::TRVector(const TRMatrix&)", "index i %d out of bounds (size: %d, this: %p)",
            I, S.NI(), this); I = 1;
  }

  fNrows = S.NJ();
  fNcols = 1;
  Set(fNrows, S.GetArray() + S.NJ() * (I - 1));
}


TRVector::TRVector(const TRVector &A, ETRMatrixCreatorsOp kop, const TRMatrix &B)
{
  Int_t NI, NJ, NK;

  switch (kop) {
  case kAxB:
    NI = A.GetNrows(); fNcols = NI;
    NJ = A.GetNcols();
    assert (NJ == B.GetNrows());
    NK = B.GetNcols(); fNrows = NK;
    Set(NI * NK);
    TCL::mxmpy(A.GetArray(), B.GetArray(), fArray, NI, NJ, NK);
    break;

  default:
    Error("TRVector(ETRMatrixCreatorsOp)", "operation %d not yet implemented", kop);
  }
}


TRVector::TRVector(Int_t nrows, const Char_t* s): TRMatrix(nrows, 1, s) {;}


TRVector::TRVector(const TRSymMatrix &S, ETRMatrixCreatorsOp kop, const TRVector &A) :
  TRMatrix(S, kop, A) {;}


TRVector::TRVector(const TRVector &A, ETRMatrixCreatorsOp kop, const TRSymMatrix &S) :
  TRMatrix(A, kop, S) {;}


TRVector::TRVector(Int_t nrows, Double_t a0, ...) : TRMatrix(nrows, 1)
{
  __VA_LIST__(a0);
}


TRVector::TRVector(const TVector3 &A) : TRMatrix(3, 1)
{
  Double_t xyz[3];
  A.GetXYZ(xyz);
  Set(3, xyz);
}


TRVector::TRVector(const StThreeVectorD &A) : TRMatrix(3, 1)
{
  Set(3, A.xyz());
}


TRVector &TRVector::operator=(const TVector3 &A)
{
  Double_t xyz[3];
  A.GetXYZ(xyz);
  Set(3, xyz);
  fNrows = 3;
  fNcols = 1;
  return *this;
}


TRVector &TRVector::operator=(const StThreeVectorD &A)
{
  Set(3, A.xyz());
  fNrows = 3;
  fNcols = 1;
  return *this;
}


ostream &operator<<(ostream &s, const TRVector &target)
{
  Int_t Nrows = target.GetNrows();
  assert(target.GetNcols() == 1);
  const Double_t* Array = target.GetArray();
  s << "Vector[" << Nrows << "] = ";// << endl;

  if (Array) for (int i = 0; i < Nrows; i++) s << Form("\t%10.3f", Array[i]);
  else s << " Empty";

  //  s << endl;
  return s;
}


void TRVector::Print(Option_t* opt) const {if (opt) {}; cout << *this << endl;}


TRVector TRVector::Cross(const TRVector &v) const
{
  assert(fNrows == 3 && fNcols == 1 && v.GetNrows() == 3 && v.GetNcols() == 1);
  TRVector out(3);
  TMath::Cross(GetArray(), v.GetArray(), out.GetArray());
  return out;
}


TRVector TRVector::Unit() const
{
  Double_t Norm = TMath::Sqrt((*this) * (*this));
  return Norm > 0 ? (*this) / Norm : *this;
}
