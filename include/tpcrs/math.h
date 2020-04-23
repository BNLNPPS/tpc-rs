#ifndef tpcrs_math_h
#define tpcrs_math_h

namespace tpcrs {

template<typename T>
int irint(T x)
{
   int i;
   if (x >= 0) {
      i = int(x + 0.5);
      if ( i & 1 && x + 0.5 == T(i) ) i--;
   } else {
      i = int(x - 0.5);
      if ( i & 1 && x - 0.5 == T(i) ) i++;
   }
   return i;
}

}

#endif
