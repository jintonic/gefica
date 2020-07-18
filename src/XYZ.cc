#include "XYZ.h"
#include "Units.h"
#include "SquarePointContact.h"
using namespace GeFiCa;

void XYZ::SetupWith(Detector &detector)
{
   Grid::SetupWith(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("SquarePointContact")) {
      SquarePointContact& spc = (SquarePointContact&) detector;
      spc.CheckConfigurations();
      GetInfoFrom(spc);
   } else {
      Error("SetupWith", "%s is not expected.", type.Data());
      Error("SetupWith", "Please use SquarePointContact detector.");
      abort();
   }

   fDetector = &detector; // for GetC to use fDetector->Bias[]
}
//_____________________________________________________________________________
//
void XYZ::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx]) return; // no need to calculate on boundaries

   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double h0=dC3m[idx];
   double h5=dC3p[idx];
   /*
   double dxm=dC1m[idx];
   double dxp=dC1p[idx];
   double dym=dC2m[idx];
   double dyp=dC2p[idx];
   double dzm=dC3m[idx];
   double dzp=dC3p[idx];
   */
   double vym,vyp,vxm,vxp,vzp,vzm;
   if (idx<N1*N2)vzm=Vp[idx];
   else vzm=Vp[idx-N1*N2];
   if (idx>=N1*N2*N3-N1*N2)vzp=Vp[idx];
   else vzp=Vp[idx+N1*N2];
   if (idx%(N1*N2)>(N1*N2)-N1-1) vyp=Vp[idx];
   else vyp=Vp[idx+N1];
   if (idx%(N1*N2)<N1)vym=Vp[idx];
   else vym=Vp[idx-N1];
   if ((idx%(N1*N2))%N1==N1-1)vxp=Vp[idx];
   else vxp=Vp[idx+1];
   if ((idx%(N1*N2))%N1==0)vxm=Vp[idx];
   else vxm=Vp[idx-1];
   double tmp= (
         Src[idx]*h0*h1*h2*h3*h4*h5*(h1+h4)*(h2+h3)*(h0+h5)/2
         +(vxp*h3+vxm*h2)*h0*h1*h4*h5*(h1+h4)*(h0+h5)
         +(vyp*h4+vym*h1)*h0*h2*h3*h5*(h0+h5)*(h2+h3)
         +(vzp*h5+vzm*h0)*h1*h2*h3*h4*(h1+h4)*(h2+h3) 
         )
      /((h0+h5)*(h1+h4)*(h2+h3)*(h0*h1*h4*h5+h0*h2*h3*h5+h1*h2*h3*h4));
   /*
      double tmp=(Src[idx]
      +(vxp/dxp+vxm/dxm)/(dxp+dxm)
      +(vyp/dyp+vym/dym)/(dyp+dym)
      +(vzp/dzp+vzm/dzm)/(dzp+dzm)
      )/(
      (1/dxp+1/dxm)/(dxp+dxm)
      +(1/dyp+1/dym)/(dyp+dym)
      +(1/dzp+1/dzm)/(dzp+dzm)
      );
      */
   double oldP=Vp[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;
   Vp[idx]=tmp;
   return;
   // update Vp for impurity-only case even if the point is undepleted
   if (fDetector->Bias[0]==fDetector->Bias[1]) { Vp[idx]=tmp; return; }

   //check depletion
   double min=vxm;
   double max=vxm;
   if(min>vxp)min=vxp;
   if (min>vyp)min=vyp;
   if (min>vym)min=vym;
   if (min>vzp)min=vzp;
   if (min>vzm)min=vzm;

   //find max
   if(max<vxp)max=vxp;
   if (max<vyp)min=vyp;
   if (max<vym)max=vym;
   if (max<vzp)max=vzp;
   if (max<vzm)max=vzm;
   //if tmp is greater or smaller than max and min, set tmp to it.
   //Vp[idx]=RelaxationFactor*(tmp-Vp[idx])+Vp[idx];
   //if need calculate depleted voltage
   if(tmp<min) {
      Vp[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      Vp[idx]=max;
      fIsDepleted[idx]=false;
   } else {
      Vp[idx]=tmp;
      fIsDepleted[idx]=true;
   }
}
//______________________________________________________________________________
//
double XYZ::GetC()
{
   //FIXME:function of integration need to be update for xyz
   return -1;

   Grid::GetC(); // calculate field excluding undepleted region

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   SquarePointContact& spc = (SquarePointContact&) *fDetector;
   double dV=spc.Bias[0]-spc.Bias[1]; if(dV<0)dV=-dV;
   double SumofElectricField=0;
   for(size_t i=0;i<GetN();i++) {
      double e1=E1[i];
      double e2=E2[i];
      double dr=dC1p[i];
      double dz=dC2p[i];
      SumofElectricField+=(e1*e1+e2*e2)*C1[i]*dr*dz;
   }
   double c=SumofElectricField*2*3.14159*epsilon/dV/dV;
   Info("GetC","%.2f pF",c/pF);
   return c;
}
//_____________________________________________________________________________
//
void XYZ::GetInfoFrom(SquarePointContact &spc)
{
   // basic setups
   double dx=detector.Width/(N1-1);
   double dy=detector.Length/(N2-1);
   double dz=detector.Height/(N3-1);

   for (size_t i=0; i<N3; i++) {
      for (size_t j=0; j<N2; j++) {
         for (size_t k=0; k<N1; k++) {
            dC1p.push_back(dx); dC1m.push_back(dx);
            dC2p.push_back(dy); dC2m.push_back(dy);
            dC3p.push_back(dz); dC3m.push_back(dz);
            C1.push_back(k*dx); C2.push_back(j*dy); C3.push_back(i*dz);
            E1.push_back(0); E2.push_back(0); E3.push_back(0);
            Et.push_back(0); Vp.push_back(0);
            fIsFixed.push_back(false); fIsDepleted.push_back(false);
         }
      }
   }

   // detailed tuning
   for (size_t i=0; i<GetN(); i++) {
      Src.push_back(-detector.GetImpurity(C3[i])*Qe/epsilon); // impurity

      // outer contact
      if (C1[i]<=0+1e-5||C1[i]>=spc.Width-1e-5
            ||C2[i]<=0+1e-5||C2[i]>=spc.Length-1e-5
            ||C3[i]<=0+1e-5||C3[i]>=spc.Height-1e-5) {
         fIsDepleted[i]=true;
         fIsFixed[i]=true;
         Vp[i]=spc.Bias[0];
         continue;
      }
      //point contact
      if(C3[i]<=spc.PointContactH&&
            C1[i]>=(spc.Width-spc.PointContactW)/2&&
            C1[i]<=(spc.Width+spc.PointContactW)/2&&
            C2[i]>=(spc.Length-spc.PointContactL)/2&&
            C2[i]<=(spc.Length+spc.PointContactL)/2) {
         fIsDepleted[i]=true;
         fIsFixed[i]=true;
         Vp[i]=spc.Bias[1];
         continue;
      }
   }
}
//___________________________________________________________________________
//
void XYZ::CalculateE()
{
   for (size_t idx=0; idx<GetN(); idx++) { 
      double dxm=dC1m[idx];
      double dxp=dC1p[idx];
      double dym=dC2m[idx];
      double dyp=dC2p[idx];
      double dzm=dC3m[idx];
      double dzp=dC3p[idx];
      double vym,vyp,vxm,vxp,vzp,vzm;
      if (idx<N1*N2) vzm=Vp[idx];
      else vzm=Vp[idx-N1*N2];
      if (idx>=N1*N2*N3-N1*N2) vzp=Vp[idx];
      else vzp=Vp[idx+N1*N2];
      if (idx%(N1*N2)>(N1*N2)-N1-1) vyp=Vp[idx];
      else vyp=Vp[idx+N1];
      if (idx%(N1*N2)<N1) vym=Vp[idx];
      else vym=Vp[idx-N1];
      if ((idx%(N1*N2))%N1==N1-1) vxp=Vp[idx];
      else vxp=Vp[idx+1];
      if ((idx%(N1*N2))%N1==0) vxm=Vp[idx];
      else vxm=Vp[idx-1];
      E1[idx]=(vxp-vxm)/(dxm+dxp);
      E2[idx]=(vyp-vym)/(dym+dyp);
      E3[idx]=(vzp-vzm)/(dzm+dzp);

      Et[idx]=sqrt(E1[idx]*E1[idx]+E2[idx]*E2[idx]+E3[idx]*E3[idx]);//overall E
   }
}
