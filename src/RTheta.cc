#include "Units.h"
#include "RTheta.h"
#include "Hemispherical.h"
using namespace GeFiCa;

void RTheta::GetBoundaryConditionFrom(Detector &detector)
{
   Grid::GetBoundaryConditionFrom(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("Hemispherical")==false) {
      Error("GetBoundaryConditionFrom", "%s is not expected. "
            "Please pass in a Hemispherical detector.", type.Data());
      abort();
   }
   Hemispherical& hemi = (Hemispherical&) detector;
   hemi.CheckConfigurations();
   fDetector = &detector; // for GetC to use fDetector->Bias[]

   double dR=hemi.Height-hemi.PointContactR;
   for (size_t i=0; i<N1; i++) {
      dC1p.push_back(dR/(N1-1)); dC1m.push_back(dR/(N1-1));
      C1.push_back(hemi.PointContactR+i*dC1p[i]);
      E1.push_back(0); Et.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-hemi.GetImpurity(C1[i])*Qe/epsilon);
   }
   dC1m[0]=0; dC1p[N1-1]=0;
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[N1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (hemi.Bias[1]-hemi.Bias[0])/(N1-1);
   for (size_t i=0; i<N1; i++) Vp.push_back(hemi.Bias[0]+slope*i);
   Vp[N1-1]=hemi.Bias[1];
}
//_____________________________________________________________________________
//
void RTheta::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return; // no need to calculate on boundaries

   double vnew=Src[idx]*dC1m[idx]*dC1p[idx]/2+(1/C1[idx]*(Vp[idx+1]-Vp[idx-1])
         +Vp[idx+1]/dC1p[idx]+Vp[idx-1]/dC1m[idx])/(1/dC1m[idx]+1/dC1p[idx]);
   vnew=RelaxationFactor*(vnew-Vp[idx])+Vp[idx];

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
void RTheta::CalculateE()
{
   for (size_t i=1; i<N1-1; i++) {
      E1[i]=(C1[i+1]*Vp[i+1]-C1[i-1]*Vp[i-1])/(dC1p[i]+dC1m[i])/C1[i];
      Et[i]=E1[i];
   }
   E1[0]=(C1[1]*Vp[1]-C1[0]*Vp[0])/dC1p[0]/C1[0]; Et[0]=E1[0];
   E1[N1-1]=(C1[N1-1]*Vp[N1-1]-C1[N1-2]*Vp[N1-2])/dC1p[N1-2]/C1[N1-2];
   Et[N1-1]=E1[N1-1];
}
