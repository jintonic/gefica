#include "Units.h"
#include "Planar1D.h"
using namespace GeFiCa;

Planar1D::Planar1D(int n) : X(n), Thickness(1*cm)
{
   SetName("p1d");
   SetTitle("1D Planar detector");
} 
//_____________________________________________________________________________
//
void Planar1D::Initialize()
{
   if (Thickness<=0) {
      Warning("Initialize", "Thickness(%.1f)<=0, set it to 1*cm", Thickness);
      Thickness=1*cm;
   }
   double steplength=(UpperBound-LowerBound)/(n-1);
   SetStepLength(steplength);
   for (int i=n; i-->0;) fC1[i]=fC1[i]+LowerBound;
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (V1-V0)/(n-1);
   for (int i=0; i<n; i++) fV[i]=V0+slope*i;
   fV[n-1]=V1;
}
//_____________________________________________________________________________
//
#include  <cmath>
bool Planar1D::Analytic()
{
   bool isConstantImpurity=true;
   for(int i=0;i+1<n;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if(!isConstantImpurity) {
      Warning("Analytic","can't handle changing impurity! Return false.");
      return false;
   }
   double d=UpperBound-LowerBound;//thickness or depth of the detector
   double a=-fImpurity[n-1]*Qe/2/epsilon;
   double b=(fV[n-1]-fV[0]-a*d*d)/d;
   double c=fV[0];
   for (int i=0; i<n; i++) fV[i] = a*fC1[i]*fC1[i]+b*fC1[i]+c;

   for (int i=1; i<n-1; i++)
      fE1[i]=(fV[i+1]-fV[i-1])/(fdC1p[i]+fdC1m[i]);
   fE1[0]=(fV[1]-fV[0])/fdC1p[0];
   fE1[n-1]=(fV[n-1]-fV[n-2])/fdC1m[n-1];

   return true;
}
//_____________________________________________________________________________
//
bool Planar1D::CalculatePotential(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return X::CalculatePotential(method);
}
