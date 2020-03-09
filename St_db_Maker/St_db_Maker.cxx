#include "TDatime.h"
#include "TTable.h"

#include "St_db_Maker/St_db_Maker.h"
#include "St_db_Maker/StValiSet.h"


Int_t  St_db_Maker::GetValidity(const TTable *tb, TDatime *const val)
{
   if (!tb)                             return -1;
   TString ts("."); ts+=tb->GetName();
   TDataSet *pa = tb->GetParent();
   const StValiSet *vs =0;
   if (pa) { vs = (StValiSet*)pa->Find(ts);}
   else    { vs = (StValiSet*)tb->Find(ts);}
   if (!vs)                             return -2;
   if (val) {
     val[0] = vs->fTimeMin;
     val[1] = vs->fTimeMax;
   }
   return vs->fVers;
}

