#include <iostream>
using namespace std;

#include "TrueCoaxial1D.h"
#include "Units.h"
using namespace GeFiCa;

void TrueCoaxial1D::Initialize()
{
   if (InnerRadius>=OuterRadius) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            InnerRadius, OuterRadius);
      return;
   }
   double steplength=(OuterRadius-InnerRadius)/(n-1);
   SetStepLength(steplength);
   for(int i=n;i-->0;)fC1[i]=fC1[i]+InnerRadius;
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (V1-V0)/(n-1);
   for (int i=0; i<n; i++) fPotential[i]=V0+slope*i;
}
//_____________________________________________________________________________
//
#include  <cmath>
bool TrueCoaxial1D::Analytic()
{
   bool isConstantImpurity=true;
   for(int i=0;i+1<n;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if(!isConstantImpurity) {
      Warning("Analytic","can't handle changing impurity! Return false.");
      return false;
   }
   double density=fImpurity[0]*Qe;
   double b=(fPotential[n-1]-fPotential[0] 
         + density*(fC1[n-1]*fC1[n-1]-fC1[0]*fC1[0])/epsilon/4)
      /(log(fC1[n-1]/fC1[0]));
   double a=fPotential[0]+density*fC1[0]*fC1[0]/epsilon/4-b*log(fC1[0]);
   for (int i=0; i<n; i++) {
      fPotential[i] = a+b*log(fC1[i])-density/4/epsilon*fC1[i]*fC1[i];

      fE1[i]=(fPotential[i+1]-fPotential[i-1])
         /(fDistanceToNext[i]+fDistanceToPrevious[i]);
   }
   return true;
}
//_____________________________________________________________________________
//
bool TrueCoaxial1D::CalculateField(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return X::CalculateField(method);
}
