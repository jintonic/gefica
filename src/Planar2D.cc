#include "Planar2D.h"
using namespace GeFiCa;

void Planar2D::Initialize()
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
   for (int i=n; i-->0;) fPotential[i]=V0;
   for(int i=0;i<n2;i++) {
      fIsFixed[i*n1]=true;
      fIsFixed[(i+1)*n1-1]=true;
      fPotential[i*n1]=V0;
      fPotential[(i+1)*n1-1]=V1;
   }
}
//_____________________________________________________________________________
//
bool Planar2D::CalculatePotential(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return XY::CalculatePotential(method);
}
