#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

Planar::Planar(const char *name, const char *title) : Detector(name, title)
{ Height=1*cm; Width=2*cm; Bias.push_back(0); Bias.push_back(1*kV); }
//_____________________________________________________________________________
//
void Planar::FillGridWithAnalyticResult()
{
   if (dC1p[0]==0) Initialize(); // setup and initialize grid if it's not done

   bool isConstantImpurity=true;
   for (int i=0;i+1<grid.N1;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if (isConstantImpurity==false) {
      Warning("FillGridWithAnalyticResult",
            "can't handle changing impurity! Abort.");
      abort();
   }
   double d=Height; //Height or depth of the detector
   double a=-fImpurity[grid.N1-1]*Qe/2/epsilon;
   double b=(V[grid.N1-1]-V[0]-a*d*d)/d;
   double c=V[0];
   for (int i=0; i<grid.N1; i++) V[i] = a*C1[i]*C1[i]+b*C1[i]+c;

   for (int i=1; i<grid.N1-1; i++)
      E1[i]=(V[i+1]-V[i-1])/(dC1p[i]+dC1m[i]);
   E1[0]=(V[1]-V[0])/dC1p[0];
   E1[grid.N1-1]=(V[grid.N1-1]-V[grid.N1-2])/dC1m[grid.N1-1];
}
//_____________________________________________________________________________
//
void Planar::ConfigureX(X& grid)
{
   if (Height<=0) {
      Error("ConfigureX", "Height(%.1fcm)<=0, abort!", Height/cm);
      abort();
   }
   size_t n1 = grid.N1;
   for (size_t i=0; i<n1; i++) {
      grid.dC1p.push_back(Height/(n1-1)); grid.dC1m.push_back(Height/(n1-1));
      grid.C1.push_back(i*dC1p[i]); grid.fIsFixed.push_back(false);
   }
   // fix 1st and last points
   grid.fIsFixed[0]=true; grid.fIsFixed[n1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (Bias[1]-Bias[0])/(n1-1);
   for (size_t i=0; i<n1-1; i++) grid.V[i]=Bias[0]+slope*i;
   grid.V[n1-1]=Bias[1];
}
//_____________________________________________________________________________
//
void Planar::ConfigureXY(XY& grid)
{
   size_t n1 = grid.N1;
   for (size_t i=0; i<n1; i++) {
      grid.dC1p.push_back(Height/(n1-1)); grid.dC1m.push_back(Height/(n1-1));
      grid.C1.push_back(i*dC1p[i]); grid.fIsFixed.push_back(false);
   }
   // fix 1st and last lines
   //size_t n1=grid.N1, n2=grid.N2, n=n1*n2;
   //for (size_t i=0;i<n;i++) {
   //   if(i>n1-1)C2[i]=C2[i-n1]+steplength2;
   //   else C2[i]=0;
   //   if(i%n1==0)C1[i]=0;
   //   else C1[i]=C1[i-1]+steplength1;

   //   E2[i]=0;
   //   dC2m[i]=steplength2;
   //   dC2p[i]=steplength2;
   //}
}
