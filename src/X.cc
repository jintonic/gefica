#include "X.h"
#include "Units.h"
using namespace GeFiCa;

void X::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return;

   double tmp = -Src[idx]*dC1m[idx]*dC1p[idx]/2 +
      (dC1p[idx]*V[idx-1]+dC1m[idx]*V[idx+1])/(dC1m[idx]+dC1p[idx]);

   tmp=RelaxationFactor*(tmp-V[idx])+V[idx];

   double min=V[idx-1], max=V[idx-1];
   if (min>V[idx+1]) min=V[idx+1];
   if (max<V[idx+1]) max=V[idx+1];

   if (tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if (tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||V.begin()==V.end()) V[idx]=tmp;
}
//_____________________________________________________________________________
//
void Planar::FillGridWithAnalyticResult(X& grid)
{
   if (TopImpurity!=BottomImpurity) {
      Error("FillGridWithAnalyticResult",
            "can't handle changing impurity! Abort.");
      abort();
   }
   double d=Height; //Height or depth of the detector
   double a=-grid.Src[grid.N1-1]/2;
   double b=(grid.V[grid.N1-1]-grid.V[0]-a*d*d)/d;
   double c=grid.V[0];
   for (size_t i=0; i<grid.N1; i++)
      grid.V[i] = a*grid.C1[i]*grid.C1[i]+b*grid.C1[i]+c;

   for (size_t i=1; i<grid.N1-1; i++)
      grid.E1[i]=(grid.V[i+1]-grid.V[i-1])/(grid.dC1p[i]+grid.dC1m[i]);
   grid.E1[0]=(grid.V[1]-grid.V[0])/grid.dC1p[0];
   grid.E1[grid.N1-1]=
      (grid.V[grid.N1-1]-grid.V[grid.N1-2])/grid.dC1m[grid.N1-1];
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
      grid.C1.push_back(i*grid.dC1p[i]); grid.fIsFixed.push_back(false);
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
