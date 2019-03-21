#include "Units.h"
#include "TrueCoaxial1D.h"
using namespace GeFiCa;

TrueCoaxial1D::TrueCoaxial1D(int n, const char *name, const char *title)
   : Rho(n, name, title), OuterR(3*cm), InnerR(0.5*cm) {};
//_____________________________________________________________________________
//
void TrueCoaxial1D::InitializeGrid()
{
   if (InnerR>=OuterR) Fatal("InitializeGrid",
         "Inner R (%.1f) >= Outer R (%.1f)! Abort!", InnerR, OuterR);
 
   double stepLength=(OuterR-InnerR)/(fN-1);
   SetStepLength(stepLength);
   for (int i=fN;i-->0;) fC1[i]=fC1[i]+InnerR;
   fIsFixed[0]=true; fIsFixed[fN-1]=true;
   double slope = (V1-V0)/(fN-1);
   for (int i=0; i<fN; i++) fV[i]=V0+slope*i;
}
//_____________________________________________________________________________
//
#include  <cmath>
using namespace std;
void TrueCoaxial1D::FillGridWithAnalyticResult()
{
   if (fdC1p[0]==0) Initialize(); // setup and initialize grid if it's not done

   bool isConstantImpurity=true;
   for (int i=0;i+1<fN;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if (!isConstantImpurity) {
      Warning("FillGridWithAnalyticResult",
            "can't handle changing impurity! Abort.");
      abort();
   }
   double density=fImpurity[0]*Qe;
   double b=(fV[fN-1]-fV[0] 
         + density*(fC1[fN-1]*fC1[fN-1]-fC1[0]*fC1[0])/epsilon/4)
      /(log(fC1[fN-1]/fC1[0]));
   double a=fV[0]+density*fC1[0]*fC1[0]/epsilon/4-b*log(fC1[0]);
   for (int i=0; i<fN; i++) {
      fV[i] = a+b*log(fC1[i])-density/4/epsilon*fC1[i]*fC1[i];
      fE1[i]=(fV[i+1]-fV[i-1])/(fdC1p[i]+fdC1m[i]);
   }
}
