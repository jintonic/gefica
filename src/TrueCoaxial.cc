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
   for (int i=fN;i-->0;) C1[i]=C1[i]+InnerR;
   fIsFixed[0]=true; fIsFixed[fN-1]=true;
   double slope = (Bias[1]-Bias[0])/(fN-1);
   for (int i=0; i<fN; i++) V[i]=Bias[0]+slope*i;
}
//_____________________________________________________________________________
//
#include  <cmath>
using namespace std;
void TrueCoaxial1D::FillGridWithAnalyticResult()
{
   if (dC1p[0]==0) Initialize(); // setup and initialize grid if it's not done

   bool isConstantImpurity=true;
   for (int i=0;i+1<fN;i++)
      if (fImpurity[i]!=fImpurity[i+1]) isConstantImpurity=false;
   if (!isConstantImpurity) {
      Warning("FillGridWithAnalyticResult",
            "can't handle changing impurity! Abort.");
      abort();
   }
   double density=fImpurity[0]*Qe;
   double b=(V[fN-1]-V[0] 
         + density*(C1[fN-1]*C1[fN-1]-C1[0]*C1[0])/epsilon/4)
      /(log(C1[fN-1]/C1[0]));
   double a=V[0]+density*C1[0]*C1[0]/epsilon/4-b*log(C1[0]);
   for (int i=0; i<fN; i++) {
      V[i] = a+b*log(C1[i])-density/4/epsilon*C1[i]*C1[i];
      E1[i]=(V[i+1]-V[i-1])/(dC1p[i]+dC1m[i]);
   }
}
