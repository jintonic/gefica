#include "TrueCoaxial2D.h"
using namespace GeFiCa;

void TrueCoaxial2D::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(OuterRadius-InnerRadius)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (InnerRadius>=OuterRadius) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            InnerRadius, OuterRadius);
      return;
   }
   double steplength=(OuterRadius-InnerRadius)/(n1-1);
   SetStepLength(steplength,3.14159265/n2);

   for(int i=n;i-->0;) fC1[i]=fC1[i]+InnerRadius;
   for (int i=0;i<n;i=i+n1)
   {
      fV[i]=V1;
      fV[i+n1-1]=V0;
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
   }
}
//_____________________________________________________________________________
//
bool TrueCoaxial2D::CalculatePotential(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return RhoPhi::CalculatePotential(method);
}
