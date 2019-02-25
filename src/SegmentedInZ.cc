#include "SegmentedInZ.h"
using namespace GeFiCa;

void SegmentedInZ::Initialize()
{
   if (InnerRadius>=OuterRadius) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            InnerRadius, OuterRadius);
      return;
   }
   double steplength=(OuterRadius-InnerRadius)/(n1-1);
   SetStepLength(steplength,2*3.14159265/n2);

   for(int i=n;i-->0;) fC1[i]=fC1[i]+InnerRadius;
   for (int i=0;i<n;i=i+n1) {
      fV[i]=V1;
      fV[i+n1-1]=V0;
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
   }
}
