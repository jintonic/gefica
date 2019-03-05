#include "Units.h"
#include "Planar1D.h"
using namespace GeFiCa;

Planar1D::Planar1D(int n, const char *name, const char *title)
   : X(n, name, title), Thickness(1*cm) {};
//_____________________________________________________________________________
//
void Planar1D::InitializeGrid()
{
   if (Thickness<=0) {
      Warning("InitializeGrid", "Thickness(%.1f)<=0, set it to 1*cm", Thickness);
      Thickness=1*cm;
   }
   // initialize fC1, fdC1p, fdC1m, fIsFixed
   SetStepLength(Thickness/(fN-1));
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[fN-1]=true;
   // linear interpolation between V0 and V1
   double slope = (V1-V0)/(fN-1);
   for (int i=0; i<fN; i++) fV[i]=V0+slope*i;
   fV[fN-1]=V1;
}
//_____________________________________________________________________________
//
#include  <cmath>
bool Planar1D::Analytic()
{
   bool isConstantImpurity=true;
   for(int i=0;i+1<fN;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if (isConstantImpurity==false) {
      Warning("Analytic","can't handle changing impurity! Return false.");
      return false;
   }
   double d=Thickness; //thickness or depth of the detector
   double a=-fImpurity[fN-1]*Qe/2/epsilon;
   double b=(fV[fN-1]-fV[0]-a*d*d)/d;
   double c=fV[0];
   for (int i=0; i<fN; i++) fV[i] = a*fC1[i]*fC1[i]+b*fC1[i]+c;

   for (int i=1; i<fN-1; i++)
      fE1[i]=(fV[i+1]-fV[i-1])/(fdC1p[i]+fdC1m[i]);
   fE1[0]=(fV[1]-fV[0])/fdC1p[0];
   fE1[fN-1]=(fV[fN-1]-fV[fN-2])/fdC1m[fN-1];

   return true;
}
