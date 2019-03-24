#include "X.h"
#include "Units.h"
using namespace GeFiCa;

void X::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return;

   double tmp = -Src[idx]*dC1m[idx]*dC1p[idx]/2 +
      (dC1p[idx]*V[idx-1]+dC1m[idx]*V[idx+1])/(dC1m[idx]+dC1p[idx]);

   tmp=RelaxationFactor*(tmp-V[idx])+V[idx];

   double min=V[idx-1], max=V[idx-1];
   if (min>V[idx+1]) min=V[idx+1];
   if (max<V[idx+1]) max=V[idx+1];

   if (tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if (tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||V.begin()==V.end()) V[idx]=tmp;
}
