#include "Units.h"
#include "RhoPhi.h"
#include "Segmented.h"
using namespace GeFiCa;

void RhoPhi::SetupWith(Detector &detector)
{
   Grid::SetupWith(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("Segmented")==false) {
      Error("SetupWith", "%s is not expected. "
            "Please pass in a Segmented detector.", type.Data());
      abort();
   }
   Segmented& sip = (Segmented&) detector;
   sip.CheckConfigurations();
   fDetector = &detector; // for GetC to use fDetector->Bias[]

   double dR=sip.Radius-sip.BoreR;
   for (size_t i=0; i<N1*N2; i++) {
      dC1p.push_back(dR/(N1-1)); dC1m.push_back(dR/(N1-1));
      dC2p.push_back(2*Pi/N2); dC2m.push_back(2*Pi/N2);
      C1.push_back(sip.BoreR+i%N1*dC1p[i]); C2.push_back(i/N1*2*Pi/N2);
      E1.push_back(0); E2.push_back(0); Et.push_back(0); Vp.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-sip.TopImpurity*Qe/epsilon);
      // fine tuning
      if (i%N1==0) { // inner surface
         dC1m[i]=0; fIsFixed[i]=true; Vp[i]=sip.Bias[0];
      } else if ((i+1)%N1==0) { // outer surface
         dC1p[i]=0; fIsFixed[i]=true;
         if (C2[i]>=2*Pi/sip.Nphi*(sip.SegmentId-1)
               && C2[i]<=2*Pi/sip.Nphi*sip.SegmentId)
            Vp[i]=sip.Bias[1]; // weighting potential for the selected segment
         else Vp[i]=sip.Bias[0]; // weighting potential for other segments
      } else
         Vp[i]=(sip.Bias[0]+sip.Bias[1])/2;
   }
}
//______________________________________________________________________________
//
void RhoPhi::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return; // no need to calculate on boundaries

   double vphip=Vp[idx+N1],vphim=Vp[idx-N1]; // v_phi_+/-
   double vrhop=Vp[idx+1],vrhom=Vp[idx-1]; // v_rho_+/-
   if (idx<N1) vphim=Vp[idx+GetN()-N1]; // phi==0 line
   if (idx>=GetN()-N1) vphip=Vp[idx-GetN()+N1]; // phi==2pi line
   double vnew = (vrhop/(dC1p[idx]*(dC1m[idx]+dC1p[idx]))
         +vrhom/(dC1m[idx]*(dC1m[idx]+dC1p[idx]))
         +vphip/C1[idx]/C1[idx]/dC2m[idx]/(dC2p[idx]+dC2m[idx])
         +vphim/C1[idx]/C1[idx]/dC2p[idx]/(dC2p[idx]+dC2m[idx])
         -Src[idx]/2+(vrhop-vrhom)/2/C1[idx]/(dC1m[idx]+dC1p[idx]))
      /(1/dC1p[idx]/(dC1m[idx]+dC1p[idx])+1/dC1m[idx]/(dC1m[idx]+dC1p[idx])
            +1/C1[idx]/C1[idx]/dC2p[idx]/(dC2p[idx]+dC2m[idx])
            +1/C1[idx]/C1[idx]/dC2m[idx]/(dC2p[idx]+dC2m[idx]));
   vnew=RelaxationFactor*(vnew-Vp[idx])+Vp[idx]; // over relex

   double vmin=vrhom; // minimal Vp around point[idx]
   if (vmin>vrhop) vmin=vrhop;
   if (vmin>vphip) vmin=vphip;
   if (vmin>vphim) vmin=vphim;

   double vmax=vrhom; // maximal Vp around point[idx]
   if (vmax<vrhop) vmax=vrhop;
   if (vmax<vphip) vmax=vphip;
   if (vmax<vphim) vmax=vphim;

   if (vnew<vmin) {
      Vp[idx]=vmin; fIsDepleted[idx]=false;
   } else if (vnew>vmax) {
      Vp[idx]=vmax; fIsDepleted[idx]=false;
   } else {
      Vp[idx]=vnew; fIsDepleted[idx]=true;
   }
   // update Vp for impurity-only case even if the point is undepleted
   if (fDetector->Bias[0]==fDetector->Bias[1]) Vp[idx]=vnew;
}
//______________________________________________________________________________
//
void RhoPhi::CalculateE()
{
   Grid::CalculateE(); // deal with E1
   for (size_t i=0; i<GetN(); i++) { // deal with E2
      E2[i]=(Vp[i+N1]-Vp[i-N1])/(dC2p[i]+dC2m[i])/C1[i];
      if (i<N1) E2[i]=(Vp[i+N1]-Vp[i])/dC2p[i]/C1[i]; // lower boundary
      if (i>GetN()-N1) E2[i]=(Vp[i]-Vp[i-N1])/dC2m[i]/C1[i]; // upper boundary
      if (i%N1==0) E1[i]=(Vp[i+1]-Vp[i])/dC1p[i]/C1[i]; // left boundary
      if ((i+1)%N1==0) E1[i]=(Vp[i]-Vp[i-1])/dC1m[i]/C1[i]; // right boundary
      Et[i]=sqrt(E1[i]*E1[i]+E2[i]*E2[i]);
   }
}
