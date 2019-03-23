#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

Planar::Planar(int n, const char *name, const char *title)
   : X(n, name, title), Thickness(1*cm) {};
//_____________________________________________________________________________
//
void Planar::InitializeGrid()
{
   if (Thickness<=0) {
      Warning("InitializeGrid", "Thickness(%.1f)<=0, set it to 1*cm", Thickness);
      Thickness=1*cm;
   }
   // initialize C1, dC1p, dC1m, fIsFixed
   SetStepLength(Thickness/(fN-1));
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[fN-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (Bias[1]-Bias[0])/(fN-1);
   for (int i=0; i<fN; i++) V[i]=Bias[0]+slope*i;
   V[fN-1]=Bias[1];
}
//_____________________________________________________________________________
//
void Planar::FillGridWithAnalyticResult()
{
   if (dC1p[0]==0) Initialize(); // setup and initialize grid if it's not done

   bool isConstantImpurity=true;
   for (int i=0;i+1<fN;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if (isConstantImpurity==false) {
      Warning("FillGridWithAnalyticResult",
            "can't handle changing impurity! Abort.");
      abort();
   }
   double d=Thickness; //thickness or depth of the detector
   double a=-fImpurity[fN-1]*Qe/2/epsilon;
   double b=(V[fN-1]-V[0]-a*d*d)/d;
   double c=V[0];
   for (int i=0; i<fN; i++) V[i] = a*C1[i]*C1[i]+b*C1[i]+c;

   for (int i=1; i<fN-1; i++)
      E1[i]=(V[i+1]-V[i-1])/(dC1p[i]+dC1m[i]);
   E1[0]=(V[1]-V[0])/dC1p[0];
   E1[fN-1]=(V[fN-1]-V[fN-2])/dC1m[fN-1];
}
