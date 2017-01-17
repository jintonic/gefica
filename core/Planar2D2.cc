#include "Planar2D2.h"
using namespace GeFiCa;

void Planar2D2::Initialize()
{
   if (XLowerBound>=XUpperBound||YLowerBound>=YUpperBound) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            XLowerBound, XUpperBound);
      return;
   }
   double steplength1=(XUpperBound-XLowerBound)/(n1-1);
   double steplength2=(YUpperBound-YLowerBound)/(n2-1);
   XY::SetStepLength(steplength1,steplength2);
   for (int i=n; i-->0;) fC1[i]=fC1[i]+XLowerBound;
   for (int i=n; i-->0;) fC2[i]=fC2[i]+YLowerBound;
   for(int i=0;i<n;i++) {
     int x=i%n1;
     int y=i/n1;
     if(x-y==-n2/2)
     {
       fIsFixed[i]=true;
       fPotential[i]=V0;
     }
     if(x-y==n2/2)
     {
       fIsFixed[i]=true;
       fPotential[i]=V1;
     }
     /*if(x-y>-n2/2)
       fIsFixed[i]=true;
     if(x-y<n2/2)
       fIsFixed[i]=true;*/
   }
}
//_____________________________________________________________________________
//
bool Planar2D2::CalculateField(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return XY::CalculateField(method);
}
