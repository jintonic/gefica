#include "Units.h"
#include "RhoPhi.h"
#include "Segmented.h"
using namespace GeFiCa;

void RhoPhi::GetBoundaryConditionFrom(Detector &detector)
{
   Grid::GetBoundaryConditionFrom(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("Segmented")==false) {
      Error("GetBoundaryConditionFrom", "%s is not expected. "
            "Please pass in a Segmented detector.", type.Data());
      abort();
   }
   Segmented& sd = (Segmented&) detector;
   sd.CheckConfigurations();
   fDetector = &detector; // for GetC to use fDetector->Bias[]

   double dR=sd.Radius-sd.BoreR;
   for (size_t i=0; i<N1*N2; i++) {
      dC1p.push_back(dR/(N1-1)); dC1m.push_back(dR/(N1-1));
      dC2p.push_back(2*Pi/N2); dC2m.push_back(2*Pi/N2);
      if (i%N1==0) dC1m[i]=0; if ((i+1)%N1==0) dC1p[i]=0;
      C1.push_back(sd.BoreR+i*dC1p[i]);
      C2.push_back(i%N1*2*Pi/N2);
      E1.push_back(0); E2.push_back(0); Et.push_back(0); Vp.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-sd.GetImpurity(C2[i])*Qe/epsilon);
   }
   dC1m[0]=0; dC1p[N1-1]=0;
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[N1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (sd.Bias[1]-sd.Bias[0])/(N1-1);
   for (size_t i=0; i<N1; i++) Vp.push_back(sd.Bias[0]+slope*i);
   Vp[N1-1]=sd.Bias[1];
   double steplength1=(Radius-BoreR)/(N1-1);
   double steplength2=2*TMath::Pi()/(N2);
   double SegmentUpperBound=2*TMath::Pi()*SegmentId/Nphi;
   double SegmentLowerBound=2*TMath::Pi()*(SegmentId-1)/Nphi;
   SetStepLength(steplength1,steplength2);
   for (size_t i=0; i<N1*N2; i++) {
      C1[i]+=BoreR;
      if(i%N1==0) {
         fIsFixed[i]=true;
         if(SegmentId==0)V[i]=Bias[1];
         else V[i]=Bias[0];
      } else if(i%N1==N1-1) { //need shift boundary as Segmented
         fIsFixed[i]=true;
         if(SegmentId==0)V[i]=Bias[0];
         else if(C2[i]<=SegmentUpperBound&&C2[i]>=SegmentLowerBound) {
            V[i]=Bias[1];
         } else if(C2[i]>=SegmentUpperBound||C2[i]<=SegmentLowerBound) {
            V[i]=Bias[0];
         }
      } else {
         fIsFixed[i]=false;
         V[i]=0;
      }
   }
}
//______________________________________________________________________________
//
void RhoPhi::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double r=C1[idx];
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double pphim,pphip,prhom,prhop;
   if(idx>=N1)pphim=V[idx-N1];
   else pphim=V[idx+GetN()-N1];
   if(idx>=GetN()-N1)pphip=V[idx-GetN()+N1];
   else pphip=V[idx+N1];
   if(idx%N1==0)prhom=V[idx];
   else prhom=V[idx-1];
   if(idx%N1==N1-1)prhop=V[idx];
   else prhop=V[idx+1];
   double tmp = (prhop/(h3*(h2+h3))+prhom/(h2*(h2+h3))
         +pphip/r/r/h4/(h1+h4)+pphim/r/r/h1/(h1+h4)
         -density/epsilon/2+(prhop-prhom)/2/r/(h2+h3))
      /(1/h3/(h2+h3)+1/h2/(h2+h3)
            +1/r/r/h1/(h1+h4)+1/r/r/h4/(h1+h4));
   double min=prhom;
   double max=prhom;
   if(min>prhop)min=prhop;
   if (min>pphip)min=pphip;
   if (min>pphim)min=pphim;

   //find max
   if(max<prhop)max=prhop;
   if (max<pphip)max=pphip;
   if (max<pphim)max=pphim;
   //if tmp is greater or smaller than max and min, set tmp to it.
   //V[idx]=RelaxationFactor*(tmp-V[idx])+V[idx];
   //if need calculate depleted voltage
   double oldP=V[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;
   if(tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||Bias[0]==Bias[1]) V[idx]=tmp;
}
//_____________________________________________________________________________
//
double RhoPhi::GetData(double x, double y, double z, double *data)
{
   //0:Impurity 1:Potential 2:E1 3:E2
   size_t idx=FindIdx(x,y,0,GetN());
   double ab=(x-C1[idx])/dC1p[idx];
   double aa=1-ab;
   double ba=(y-C2[idx])/dC2p[idx];
   double bb=1-ba;
   double tar0,tar1,tar2,tar3;
   tar3=-1;
   tar0=data[idx];
   if((idx%N1)==N1-1) {
      tar1=data[idx+1-N1];
      tar3=data[idx+1];
   } else {
      tar1=data[idx+1];
   }
   if(idx>GetN()-N1){tar2=0;tar3=0;}
   else {tar2=data[idx+N1];}
   if (tar3==-1)tar3=data[idx+N1+1];
   return (tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb;
}
//______________________________________________________________________________
//
