#ifndef St_SurveyC_h
#define St_SurveyC_h

#include "tpcrs/config_structs.h"
#include "Survey.h"
#include "TGeoMatrix.h"
struct St_SurveyC : tpcrs::IConfigStruct {
  virtual Survey_st* Struct(int i = 0) const = 0;

  void Initialize()
  {
    UInt_t N = GetNRows();
    fRotations = new TGeoHMatrix*[N];
    for (UInt_t i = 0; i < N; i++) {
      fRotations[i] = new TGeoHMatrix;
      TGeoHMatrix &rot = *fRotations[i];
      if (N == 1) rot.SetName(GetName().c_str());
      else        rot.SetName(Form("%s_%i",GetName().c_str(),i+1));
      rot.SetRotation(Rotation(i));
      rot.SetTranslation(Translation(i));
      Normalize(rot);
      assert(TMath::Abs(rot.Determinant())-1 < 1.e-3);
    }
  }

  virtual  ~St_SurveyC();
  Int_t 	Id(Int_t i = 0) 	const {return Struct(i)->Id;}
  Double_t 	r00(Int_t i = 0) 	const {return Struct(i)->r00;} // 0
  Double_t 	r01(Int_t i = 0) 	const {return Struct(i)->r01;} // 1
  Double_t 	r02(Int_t i = 0) 	const {return Struct(i)->r02;} // 2
  Double_t 	r10(Int_t i = 0) 	const {return Struct(i)->r10;} // 3
  Double_t 	r11(Int_t i = 0) 	const {return Struct(i)->r11;} // 4
  Double_t 	r12(Int_t i = 0) 	const {return Struct(i)->r12;} // 5
  Double_t 	r20(Int_t i = 0) 	const {return Struct(i)->r20;} // 6
  Double_t 	r21(Int_t i = 0) 	const {return Struct(i)->r21;} // 7
  Double_t 	r22(Int_t i = 0) 	const {return Struct(i)->r22;} // 8
  Double_t 	t0(Int_t i = 0) 	const {return Struct(i)->t0;}
  Double_t 	t1(Int_t i = 0) 	const {return Struct(i)->t1;}
  Double_t 	t2(Int_t i = 0) 	const {return Struct(i)->t2;}
  Double_t 	sigmaRotX(Int_t i = 0) 	const {return Struct(i)->sigmaRotX;}
  Double_t 	sigmaRotY(Int_t i = 0) 	const {return Struct(i)->sigmaRotY;}
  Double_t 	sigmaRotZ(Int_t i = 0) 	const {return Struct(i)->sigmaRotZ;}
  Double_t 	sigmaTrX(Int_t i = 0) 	const {return Struct(i)->sigmaTrX;}
  Double_t 	sigmaTrY(Int_t i = 0) 	const {return Struct(i)->sigmaTrY;}
  Double_t 	sigmaTrZ(Int_t i = 0) 	const {return Struct(i)->sigmaTrZ;}
  Char_t* 	comment(Int_t i = 0) 	const {return Struct(i)->comment;}
  void          GetAngles(Double_t &phi, Double_t &the, Double_t &psi, Int_t i = 0);
  const Double_t  *Rotation(Int_t i = 0)     const {return &Struct(i)->r00;} 
  const Double_t  *Translation(Int_t i = 0)  const {return &Struct(i)->t0;} 
  const TGeoHMatrix  &GetMatrix(Int_t i = 0);
  const TGeoHMatrix  &GetMatrix4Id(Int_t id);
  const TGeoHMatrix  &GetMatrixR(Int_t i); // ignoring rotation alpha and beta
  const Double_t *r(Int_t i = 0)        const {return &Struct(i)->r00;}
  const Double_t *t(Int_t i = 0)        const {return &Struct(i)->t0;}
  static void Normalize(TGeoHMatrix &rot);
  static Double_t IsOrtogonal(const Double_t *r);
 protected:
  St_SurveyC();
 private:
  TGeoHMatrix  **fRotations;
};
#endif
