#include "Units.h"
#include "TrueCoaxial1D.h"
using namespace GeFiCa;

TrueCoaxial1D::TrueCoaxial1D(int n, const char *name, const char *title)
   : Rho(n, name, title), OuterRadius(3*cm), InnerRadius(0.5*cm) {};
//_____________________________________________________________________________
//
void TrueCoaxial1D::Initialize()
{
   if (InnerRadius>=OuterRadius) Fatal("Initialize",
         "Inner R (%.1f) >= Outer R (%.1f)! Abort!", InnerRadius, OuterRadius);
 
   double stepLength=(OuterRadius-InnerRadius)/(n-1);
   SetStepLength(stepLength);
   for (int i=n;i-->0;) fC1[i]=fC1[i]+InnerRadius;
   fIsFixed[0]=true; fIsFixed[n-1]=true;
   double slope = (V1-V0)/(n-1);
   for (int i=0; i<n; i++) fV[i]=V0+slope*i;
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
   double b=(fV[n-1]-fV[0] 
         + density*(fC1[n-1]*fC1[n-1]-fC1[0]*fC1[0])/epsilon/4)
      /(log(fC1[n-1]/fC1[0]));
   double a=fV[0]+density*fC1[0]*fC1[0]/epsilon/4-b*log(fC1[0]);
   for (int i=0; i<n; i++) {
      fV[i] = a+b*log(fC1[i])-density/4/epsilon*fC1[i]*fC1[i];

      fE1[i]=(fV[i+1]-fV[i-1])
         /(fdC1p[i]+fdC1m[i]);
   }
   return true;
}
