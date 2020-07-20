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

   double dxp=dC1p[idx]!=0?dC1p[idx]:dC1m[idx];
   double dxm=dC1m[idx]!=0?dC1m[idx]:dC1p[idx];
   double dyp=dC2p[idx]!=0?dC2p[idx]:dC2m[idx];
   double dym=dC2m[idx]!=0?dC2m[idx]:dC2p[idx];
   double dzp=dC3p[idx]!=0?dC3p[idx]:dC3m[idx];
   double dzm=dC3m[idx]!=0?dC3m[idx]:dC3p[idx];
 
   double vxp,vxm,vyp,vym,vzp,vzm;
   if (idx<N1*N2) vzm=Vp[idx]; // bottom boundary
   else vzm=Vp[idx-N1*N2];
   if (idx>=N1*N2*N3-N1*N2) vzp=Vp[idx]; // top boundary
   else vzp=Vp[idx+N1*N2];
   if (idx%(N1*N2)>(N1*N2)-N1-1) vyp=Vp[idx]; // back boundary
   else vyp=Vp[idx+N1];
   if (idx%(N1*N2)<N1) vym=Vp[idx]; // front boundary
   else vym=Vp[idx-N1];
   if ((idx%(N1*N2))%N1==N1-1) vxp=Vp[idx]; // right boundary
   else vxp=Vp[idx+1];
   if ((idx%(N1*N2))%N1==0) vxm=Vp[idx]; // left boundary
   else vxm=Vp[idx-1];

   // update potential
   double vnew = ( Src[idx]/2 + (vxp/dxp+vxm/dxm)/(dxm+dxp)
         + (vyp/dyp+vym/dym)/(dym+dyp) + (vzp/dzp+vzm/dzm)/(dzp+dzm) )
      / ( (1/dxp+1/dxm)/(dxp+dxm)
            + (1/dyp+1/dym)/(dyp+dym) + (1/dzp+1/dzm)/(dzp+dzm) );
   vnew = RelaxationFactor*(vnew-Vp[idx])+Vp[idx]; // over relax

   // update Vp for impurity-only case even if the point is undepleted
   if (fDetector->Bias[0]==fDetector->Bias[1]) { Vp[idx]=vnew; return; }

   //check depletion
   fIsDepleted[idx]=false; // default
   //find minimal potential in all neighboring points
   double vmin=vxp; // minimal Vp around point[idx]
   if (vmin>vxm)vmin=vxm;
   if (vmin>vyp)vmin=vyp;
   if (vmin>vym)vmin=vym;
   if (vmin>vzp)vmin=vzp;
   if (vmin>vzm)vmin=vzm;
   //find maximal potential in all neighboring points
   double vmax=vxp; // maximal Vp around point[idx]
   if (vmax<vxm)vmax=vxm;
   if (vmax<vyp)vmax=vyp;
   if (vmax<vym)vmax=vym;
   if (vmax<vzp)vmax=vzp;
   if (vmax<vzm)vmax=vzm;
   //if vnew is greater or smaller than max and min, set vnew to it.
   if (vnew<vmin) Vp[idx]=vmin;
   else if(vnew>vmax) Vp[idx]=vmax;
   else { Vp[idx]=vnew; fIsDepleted[idx]=true; } // vmin<vnew<vmax
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
   double dx=spc.Width/(N1-1), dy=spc.Length/(N2-1), dz=spc.Height/(N3-1);

   for (size_t i=0; i<N3; i++) {
      for (size_t j=0; j<N2; j++) {
         for (size_t k=0; k<N1; k++) {
            C1.push_back(k*dx); C2.push_back(j*dy); C3.push_back(i*dz);
            dC1p.push_back(dx); dC1m.push_back(dx);
            dC2p.push_back(dy); dC2m.push_back(dy);
            dC3p.push_back(dz); dC3m.push_back(dz);
            E1.push_back(0); E2.push_back(0); E3.push_back(0);
            Et.push_back(0); Vp.push_back(0);
            fIsFixed.push_back(false); fIsDepleted.push_back(false);
            Src.push_back(-spc.GetImpurity(C3.back())*Qe/epsilon);

            if (k==0) { // left most surface
               dC1m.back()=0; E2.back()=0; E3.back()=0; Vp.back()=spc.Bias[1];
               fIsFixed.back()=true; fIsDepleted.back()=true;
            }
            if (k==N1-1) { // right most surface
               dC1p.back()=0; E2.back()=0; E3.back()=0; Vp.back()=spc.Bias[1];
               fIsFixed.back()=true; fIsDepleted.back()=true;
            }
            if (j==0) { // front most surface
               dC2m.back()=0; E1.back()=0; E3.back()=0; Vp.back()=spc.Bias[1];
               fIsFixed.back()=true; fIsDepleted.back()=true;
            }
            if (j==N2-1) { // back most surface
               dC2p.back()=0; E1.back()=0; E3.back()=0; Vp.back()=spc.Bias[1];
               fIsFixed.back()=true; fIsDepleted.back()=true;
            }
            if (i==0) { // top surface
               dC3m.push_back(0);
               // wrap around
               if (C1.back()<=(spc.Width-spc.WrapAroundW)/2
                     || C1.back()>=(spc.Width+spc.WrapAroundW)/2) {
                  E1.back()=0; E2.back()=0; Vp.back()=spc.Bias[1];
                  fIsFixed.back()=true; fIsDepleted.back()=true;
               }
               if (C2.back()<=(spc.Length-spc.WrapAroundL)/2
                     || C2.back()>=(spc.Length+spc.WrapAroundL)/2) {
                  E1.back()=0; E2.back()=0; Vp.back()=spc.Bias[1];
                  fIsFixed.back()=true; fIsDepleted.back()=true;
               }
               // point contact
               if (C3.back()<=spc.PointContactH
                     && C1.back()>=(spc.Width-spc.PointContactW)/2
                     && C1.back()<=(spc.Width+spc.PointContactW)/2
                     && C2.back()>=(spc.Length-spc.PointContactL)/2
                     && C2.back()<=(spc.Length+spc.PointContactL)/2) {
                  fIsFixed.back()=true; fIsDepleted.back()=true;
                  Vp.back()=spc.Bias[0];
               }
            }
            if (i==N3-1) { // bottom electrode
               dC3p.push_back(0);
               E1.back()=0; E2.back()=0; Vp.back()=spc.Bias[1];
               fIsFixed.back()=true; fIsDepleted.back()=true;
            }
         }
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

      Et[idx]=sqrt(E1[idx]*E1[idx]+E2[idx]*E2[idx]+E3[idx]*E3[idx]);
   }
}
