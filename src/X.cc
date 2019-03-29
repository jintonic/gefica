#include "X.h"
#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

void X::GetBoundaryConditionFrom(Detector &detector)
{
   Grid::GetBoundaryConditionFrom(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("Planar")==false) {
      Error("GetBoundaryConditionFrom", "%s is not expected. "
            "Please pass in a Planar detector.", type.Data());
      abort();
   }
   Planar& planar = (Planar&) detector;
   planar.CheckConfigurations();

   for (size_t i=0; i<N1; i++) {
      dC1p.push_back(planar.Height/(N1-1));
      dC1m.push_back(planar.Height/(N1-1));
      C1.push_back(i*dC1p[i]);
      E1.push_back(0); Et.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-planar.GetImpurity(C1[i])*Qe/epsilon);
   }
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[N1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (planar.Bias[1]-planar.Bias[0])/(N1-1);
   for (size_t i=0; i<N1; i++) Vp.push_back(planar.Bias[0]+slope*i);
   Vp[N1-1]=planar.Bias[1];
}
//_____________________________________________________________________________
//
void X::SolveAnalytically()
{
   Grid::SolveAnalytically(); // check if impurity is constant
   double h=C1[N1-1]-C1[0];
   double a=-Src[0]/2;
   double b=(Vp[N1-1]-Vp[0]-a*h*h)/h;
   double c=Vp[0];
   for (size_t i=0; i<N1; i++) Vp[i] = a*C1[i]*C1[i]+b*C1[i]+c;
   CalculateE();
}
//_____________________________________________________________________________
//
double X::GetC()
{
   Info("GetC","Start...");
   SuccessiveOverRelax(); // identify undepleted region
   std::vector<double>original=Src; // save original rho/epsilon
   for (size_t i=0; i<GetN(); i++) {
      Src[i]=0; // set impurity to zero
      if (fIsDepleted[i]==false) fIsFixed[i]=true; // fix undepleted points
   }
   SuccessiveOverRelax(); // calculate potential without impurity
   Src=original; // set impurity back

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   double dV = Vp[N1-1]-Vp[0]; if (dV<0) dV=-dV;
   double integral=0;
   for (size_t i=0; i<GetN(); i++) {
      integral+=E1[i]*E1[i]*dC1p[i];
      if (!fIsDepleted[i]) fIsFixed[i]=false; // release undepleted points
   }
   double c=integral*epsilon/dV/dV;
   Info("GetC","%.2f pF",c/pF);
   return c;
}
//_____________________________________________________________________________
//
void X::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return; // no need to calculate on boundaries

   // calculate Vp[idx] from Vp[idx-1] and Vp[idx+1]
   double vnew = Src[idx]*dC1m[idx]*dC1p[idx]/2 +
      (dC1p[idx]*Vp[idx-1]+dC1m[idx]*Vp[idx+1])/(dC1m[idx]+dC1p[idx]);
   vnew = RelaxationFactor*(vnew-Vp[idx]) + Vp[idx]; // over relax

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
//_____________________________________________________________________________
//
void X::CalculateE()
{
   for (size_t i=1; i<N1-1; i++) {
      E1[i]=(Vp[i+1]-Vp[i-1])/(dC1p[i]+dC1m[i]); Et[i]=E1[i];
   }
   E1[0]=(Vp[1]-Vp[0])/dC1p[0]; Et[0]=E1[0];
   E1[N1-1]=(Vp[N1-1]-Vp[N1-2])/dC1m[N1-1]; Et[N1-1]=E1[N1-1];
}
