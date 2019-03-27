#include "X.h"
#include "Units.h"
#include "Detector.h"
using namespace GeFiCa;

void X::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return;

   double tmp = -Src[idx]*dC1m[idx]*dC1p[idx]/2 +
      (dC1p[idx]*Vp[idx-1]+dC1m[idx]*Vp[idx+1])/(dC1m[idx]+dC1p[idx]);

   tmp=RelaxationFactor*(tmp-Vp[idx])+Vp[idx];

   double min=Vp[idx-1], max=Vp[idx-1];
   if (min>Vp[idx+1]) min=Vp[idx+1];
   if (max<Vp[idx+1]) max=Vp[idx+1];

   if (tmp<min) {
      Vp[idx]=min;
      fIsDepleted[idx]=false;
   } else if (tmp>max) {
      Vp[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||Vp.begin()==Vp.end()) Vp[idx]=tmp;
}
//_____________________________________________________________________________
//
void X::FillGridWithAnalyticResult()
{
//   if (TopImpurity!=BottomImpurity) {
//      Error("FillGridWithAnalyticResult",
//            "can't handle changing impurity! Abort.");
//      abort();
//   }
//   double d=Height; //Height or depth of the detector
//   double a=-Src[N1-1]/2;
//   double b=(Vp[N1-1]-Vp[0]-a*d*d)/d;
//   double c=Vp[0];
//   for (size_t i=0; i<N1; i++)
//      Vp[i] = a*C1[i]*C1[i]+b*C1[i]+c;
//
//   for (size_t i=1; i<N1-1; i++)
//      E1[i]=(Vp[i+1]-Vp[i-1])/(dC1p[i]+dC1m[i]);
//   E1[0]=(Vp[1]-Vp[0])/dC1p[0];
//   E1[N1-1]=(Vp[N1-1]-Vp[N1-2])/dC1m[N1-1];
}
//_____________________________________________________________________________
//
void X::SetBoundaryCondition(Detector *detector)
{
   Printf("%s",detector->ClassName());
//   if (Height<=0) {
//      Error("ConfigureX", "Height(%.1fcm)<=0, abort!", Height/cm);
//      abort();
//   }
//   size_t n1 = N1;
//   for (size_t i=0; i<n1; i++) {
//      dC1p.push_back(Height/(n1-1)); dC1m.push_back(Height/(n1-1));
//      C1.push_back(i*dC1p[i]); fIsFixed.push_back(false);
//   }
//   // fix 1st and last points
//   fIsFixed[0]=true; fIsFixed[n1-1]=true;
//   // linear interpolation between Bias[0] and Bias[1]
//   double slope = (Bias[1]-Bias[0])/(n1-1);
//   for (size_t i=0; i<n1-1; i++) Vp[i]=Bias[0]+slope*i;
//   Vp[n1-1]=Bias[1];
}
//_____________________________________________________________________________
//
