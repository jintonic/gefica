#include "Sphere.h"
#include <iostream>
using namespace GeFiCa;
using namespace std;

void Sphere::Initialize()
{
   if (InnerRadius>=OuterRadius) {
      Warning("Initialize",
            "inner radius (%f) >= outer radius (%f)! No grid is created!",
            InnerRadius, OuterRadius);
      return;
   }
   double steplength=(OuterRadius-InnerRadius)/(n1-1);
   SetStepLength(steplength,3.14159265/n2,3.14159265*2/n3);
   for(int i=n;i-->0;) fC1[i]=fC1[i]+InnerRadius;
   for(int i=n;i-->0;) fC2[i]=fC2[i]+3.14159265/2/n2;

   for (int i=0; i<n; i+=n1) {
      fPotential[i]=V0;
      fPotential[i+n1-1]=V1;
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
   }
}
//_____________________________________________________________________________
//
bool Sphere::CalculatePotential(EMethod method)
{
   if(!fIsLoaded) Initialize();
   return X::CalculatePotential(method);
}
