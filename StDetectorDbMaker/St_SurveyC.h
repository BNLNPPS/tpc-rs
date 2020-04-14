#ifndef St_SurveyC_h
#define St_SurveyC_h

#include "tpcrs/config_structs.h"
#include "Survey.h"
#include "TGeoMatrix.h"
struct St_SurveyC : tpcrs::IConfigStruct {
  virtual Survey_st* Struct(int i = 0) const = 0;

  void Initialize()
  {
    unsigned int N = GetNRows();
    fRotations = new TGeoHMatrix*[N];
    for (unsigned int i = 0; i < N; i++) {
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
  int 	Id(int i = 0) 	const {return Struct(i)->Id;}
  double 	r00(int i = 0) 	const {return Struct(i)->r00;} // 0
  double 	r01(int i = 0) 	const {return Struct(i)->r01;} // 1
  double 	r02(int i = 0) 	const {return Struct(i)->r02;} // 2
  double 	r10(int i = 0) 	const {return Struct(i)->r10;} // 3
  double 	r11(int i = 0) 	const {return Struct(i)->r11;} // 4
  double 	r12(int i = 0) 	const {return Struct(i)->r12;} // 5
  double 	r20(int i = 0) 	const {return Struct(i)->r20;} // 6
  double 	r21(int i = 0) 	const {return Struct(i)->r21;} // 7
  double 	r22(int i = 0) 	const {return Struct(i)->r22;} // 8
  double 	t0(int i = 0) 	const {return Struct(i)->t0;}
  double 	t1(int i = 0) 	const {return Struct(i)->t1;}
  double 	t2(int i = 0) 	const {return Struct(i)->t2;}
  double 	sigmaRotX(int i = 0) 	const {return Struct(i)->sigmaRotX;}
  double 	sigmaRotY(int i = 0) 	const {return Struct(i)->sigmaRotY;}
  double 	sigmaRotZ(int i = 0) 	const {return Struct(i)->sigmaRotZ;}
  double 	sigmaTrX(int i = 0) 	const {return Struct(i)->sigmaTrX;}
  double 	sigmaTrY(int i = 0) 	const {return Struct(i)->sigmaTrY;}
  double 	sigmaTrZ(int i = 0) 	const {return Struct(i)->sigmaTrZ;}
  char* 	comment(int i = 0) 	const {return Struct(i)->comment;}
  void          GetAngles(double &phi, double &the, double &psi, int i = 0);
  const double  *Rotation(int i = 0)     const {return &Struct(i)->r00;} 
  const double  *Translation(int i = 0)  const {return &Struct(i)->t0;} 
  const TGeoHMatrix  &GetMatrix(int i = 0);
  const TGeoHMatrix  &GetMatrix4Id(int id);
  const TGeoHMatrix  &GetMatrixR(int i); // ignoring rotation alpha and beta
  const double *r(int i = 0)        const {return &Struct(i)->r00;}
  const double *t(int i = 0)        const {return &Struct(i)->t0;}
  static void Normalize(TGeoHMatrix &rot);
  static double IsOrtogonal(const double *r);
 protected:
  St_SurveyC();
 private:
  TGeoHMatrix  **fRotations;
};
#endif
