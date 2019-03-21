#include "Units.h"
#include "Sphere1D.h"
using namespace GeFiCa;

Sphere1D::Sphere1D(int n, const char *name, const char *title)
   : R(n, name, title), InnerR(0.3*cm), OuterR(3*cm) {};
//_____________________________________________________________________________
//
void Sphere1D::InitializeGrid()
{
   if (OuterR<=InnerR) Fatal("InitializeGrid",
         "Inner R (%.1f) >= outer R (%.1f)! Abort!", InnerR, OuterR);

   double stepLength=(OuterR-InnerR)/(fN-1);
   SetStepLength(stepLength);
   for (int i=fN;i-->0;) C1[i]=C1[i]+InnerR;
   fIsFixed[0]=true; fIsFixed[fN-1]=true;
   double slope = (Bias[1]-Bias[0])/(fN-1);
   for (int i=0; i<fN; i++) V[i]=Bias[0]+slope*i;
}
//_____________________________________________________________________________
//
void Sphere1D::FillGridWithAnalyticResult()
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
   double c1=(Bias[1]-Bias[0] + density/epsilon/6*(C1[fN-1]*C1[fN-1]-C1[0]*C1[0]))
      /(1/C1[fN-1]-1/C1[0]);
   double c2=Bias[0]+density/epsilon/6*C1[0]*C1[0]-c1/C1[0];
   for (int i=0; i<fN; i++) {
      V[i] = -density/6/epsilon*C1[i]*C1[i]+c1/C1[i]+c2;
      // Fixme:
      if (i!=0||i!=fN-1)E1[i]=(V[i+1]-V[i-1])/(dC1p[i]+dC1m[i]);
   }
}
