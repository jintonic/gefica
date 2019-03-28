#include "X.h"
#include "Units.h"
#include "Planar.h"
using namespace GeFiCa;

void X::GetBoundaryConditionFrom(Detector &detector)
{
   if (GetN()>0) { // this function can only be called once
      Warning("GetBoundaryConditionFrom", "has been called. Do nothing.");
      return;
   }
   TString type(detector.ClassName());
   if (type.Contains("Planar")==false) {
      Error("GetBoundaryConditionFrom", "%s is not expected. "
            "Please pass in a GeFiCa::Planar detector.", type.Data());
      abort();
   }

   if (detector.Height<=0) {
      Error("GetBoundaryConditionFrom",
            "Height(%.1fcm)<=0, abort!", detector.Height/cm);
      abort();
   }

   for (size_t i=0; i<N1; i++) {
      dC1p.push_back(detector.Height/(N1-1));
      dC1m.push_back(detector.Height/(N1-1));
      C1.push_back(i*dC1p[i]);
      E1.push_back(0); Et.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-detector.GetImpurity(C1[i])*Qe/epsilon);
   }
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[N1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (detector.Bias[1]-detector.Bias[0])/(N1-1);
   for (size_t i=0; i<N1; i++) {
      Vp.push_back(detector.Bias[0]+slope*i);
   }
   Vp[N1-1]=detector.Bias[1];

   fDetector = &detector; // save it for other uses
}
//_____________________________________________________________________________
//
void X::SolveAnalytically()
{
   if (fDetector==0) {
      Error("SolveAnalytically", "Grid is not ready. "
            "Please call GetBoundaryConditionFrom(Detector&) first.");
      abort();
   }
   if (fDetector->TopImpurity!=fDetector->BottomImpurity) {
      Error("SolveAnalytically", "can't handle changing impurity.");
      abort();
   }
   double h=fDetector->Height;
   double a=-Src[N1-1]/2;
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
   double dV = fDetector->Bias[1]-fDetector->Bias[0]; if (dV<0) dV=-dV;
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

   // over relax
   double vnew = Src[idx]*dC1m[idx]*dC1p[idx]/2 +
      (dC1p[idx]*Vp[idx-1]+dC1m[idx]*Vp[idx+1])/(dC1m[idx]+dC1p[idx]);
   vnew = RelaxationFactor*(vnew-Vp[idx]) + Vp[idx];

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
   if (fDetector->Bias[0]==fDetector->Bias[1]) Vp[idx]=vnew;
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
//_____________________________________________________________________________
//
double X::GetData(const std::vector<double> &data,
      double x, double y, double z) const
{
   //     |<---dC1m[idx]--->|
   //     +---r1---+---r2---+
   // C1[idx-1]    x      C1[idx]
   size_t idx=GetIdxOfPointToTheRightOf(x);
   if (idx==N1) return data[idx];
   double r2=(C1[idx]-x)/dC1m[idx];
   double r1=1-r2;
   double xval=data[idx]*r1+data[idx-1]*r2;
   if (gDebug>0) Info("GetData","data(x=%.1f)=%.1f C1[%zu]=%.1f, "
         "r1=%.1f, r2=%0.1f", x,xval,idx,C1[idx],r1,r2);
   return xval;
}
