#include "Rho.h"
#include "Units.h"
#include "TrueCoaxial.h"
using namespace GeFiCa;

void Rho::GetBoundaryConditionFrom(Detector &detector)
{
   Grid::GetBoundaryConditionFrom(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("TrueCoaxial")==false) {
      Error("GetBoundaryConditionFrom", "%s is not expected. "
            "Please pass in a TrueCoaxial detector.", type.Data());
      abort();
   }
   TrueCoaxial& coaxial = (TrueCoaxial&) detector;
   coaxial.CheckConfigurations();
   fDetector = &detector; // for GetC to use fDetector->Bias[]

   double dR=coaxial.Radius-coaxial.BoreR;
   for (size_t i=0; i<N1; i++) {
      dC1p.push_back(dR/(N1-1)); dC1m.push_back(dR/(N1-1));
      C1.push_back(coaxial.BoreR+i*dC1p[i]);
      E1.push_back(0); Et.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-coaxial.GetImpurity(C1[i])*Qe/epsilon);
   }
   dC1m[0]=0; dC1p[N1-1]=0;
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[N1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (coaxial.Bias[1]-coaxial.Bias[0])/(N1-1);
   for (size_t i=0; i<N1; i++) Vp.push_back(coaxial.Bias[0]+slope*i);
   Vp[N1-1]=coaxial.Bias[1];
}
//______________________________________________________________________________
//
void Rho::SolveAnalytically()
{
   Grid::SolveAnalytically(); // check if impurity is constant
   // https://www.wolframalpha.com/input/?i=(V%27(r)r)%27%2Fr%3Da
   double c1 = (Vp[N1-1]-Vp[0] + Src[0]*(C1[N1-1]*C1[N1-1]-C1[0]*C1[0])/4);
   c1/=log(C1[N1-1]/C1[0]);
   double c2 = Vp[N1-1]*log(C1[0]) - Vp[0]*log(C1[N1-1]) + 
      Src[0]*(C1[N1-1]*C1[N1-1]*log(C1[0])-C1[0]*C1[0]*log(C1[N1-1]))/4;
   c2/=log(C1[0])-log(C1[N1-1]);
   for (size_t i=0; i<N1; i++)
      Vp[i] = -Src[0]*C1[i]*C1[i]/4 + c1*log(C1[i]) + c2;
   CalculateE();
}
//______________________________________________________________________________
//
double Rho::GetC()
{
   Grid::GetC(); // calculate field excluding undepleted region

   double dV = fDetector->Bias[1]-fDetector->Bias[0]; if (dV<0) dV=-dV;
   double integral=0;
   for (size_t i=0; i<GetN(); i++) {
      integral+=E1[i]*E1[i]*dC1p[i]*C1[i];
      if (!fIsDepleted[i]) fIsFixed[i]=false; // release undepleted points
   }
   double c=integral*2*Pi*epsilon/dV/dV;
   Info("GetC","%.2f pF/cm",c/pF*cm);
   return c;
}
//______________________________________________________________________________
//
void Rho::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return; // no need to calculate on boundaries

   // calculate Vp[idx] from Vp[idx-1] and Vp[idx+1]
   double vnew = Src[idx]*dC1p[idx]*dC1m[idx]/2
         + ((Vp[idx+1]-Vp[idx-1])/C1[idx]/2
               + Vp[idx+1]/dC1p[idx]+Vp[idx-1]/dC1m[idx])
         /(1/dC1m[idx]+1/dC1p[idx]);
   vnew=RelaxationFactor*(vnew-Vp[idx])+Vp[idx]; // over relax

   // check depletion and update Vp[idx] accordingly
   double min=Vp[idx-1], max=Vp[idx-1];
   if (min>Vp[idx+1]) min=Vp[idx+1];
   if (max<Vp[idx+1]) max=Vp[idx+1];
   if (vnew<min) {
      fIsDepleted[idx]=false; Vp[idx]=min;
   } else if (vnew>max) {
      fIsDepleted[idx]=false; Vp[idx]=max;
   } else {
      fIsDepleted[idx]=true; Vp[idx]=vnew;
   }

   // update Vp for impurity-only case even if the point is undepleted
   if (Vp[0]==Vp[N1-1]) Vp[idx]=vnew;
}
