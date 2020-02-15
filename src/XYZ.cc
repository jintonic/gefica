#include "XYZ.h"
#include "Units.h"
using namespace GeFiCa;

void XYZ::SetupWith(Detector &detector)
{
   Grid::SetupWith(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("SquarePointContact")) {
      SquarePointContact& spc = (SquarePointContact&) detector;
      spc.CheckConfigurations();
      GetInfoFrom(spc);
   }
   else {
      //maybe add auto configure name list 
      Error("SetupWith", "%s is not expected.", type.Data());
      Error("SetupWith", "Please use "
            "PointContact detector.");
      abort();
   }

   fDetector = &detector; // for GetC to use fDetector->Bias[]
}
//_____________________________________________________________________________
//
void XYZ::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx])return;
   double density=Src[idx];
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double h0=dC3m[idx];
   double h5=dC3p[idx];
   double pym,pyp,pxm,pxp,pzp,pzm;
   if(idx<N1*N2)pzm=Vp[idx];
   else pzm=Vp[idx-N1*N2];
   if(idx>=N1*N2*N3-N1*N2)pzp=Vp[idx];
   else pzp=Vp[idx+N1*N2];
   if(idx%(N1*N2)>(N1*N2)-N1-1) pyp=Vp[idx];
   else pyp=Vp[idx+N1];
   if(idx%(N1*N2)<N1)pym=Vp[idx];
   else pym=Vp[idx-N1];
   if((idx%(N1*N2))%N1==N1-1)pxp=Vp[idx];
   else pxp=Vp[idx+1];
   if((idx%(N1*N2))%N1==0)pxm=Vp[idx];
   else pxm=Vp[idx-1];

   double tmp= (
         density*h0*h1*h2*h3*h4*h5*(h1+h4)*(h2+h3)*(h0+h5)/2
         +(pxp*h3+pxm*h2)*h0*h1*h4*h5*(h1+h4)*(h0+h5)
         +(pyp*h4+pym*h1)*h0*h2*h3*h5*(h0+h5)*(h2+h3)
         +(pzp*h5+pzm*h0)*h1*h2*h3*h4*(h1+h4)*(h2+h3)	
         )
      /((h0+h5)*(h1+h4)*(h2+h3)*(h0*h1*h4*h5+h0*h2*h3*h5+h1*h2*h3*h4));
   

   // update Vp for impurity-only case even if the point is undepleted
   if (fDetector->Bias[0]==fDetector->Bias[1]) { Vp[idx]=tmp; return; }


   //check depletion
   double min=pxm;
   double max=pxm;
   if(min>pxp)min=pxp;
   if (min>pyp)min=pyp;
   if (min>pym)min=pym;
   if (min>pzm)min=pzm;
   if (min>pzm)min=pzm;

   //find max
   if(max<pxp)max=pxp;
   if (max<pyp)min=pyp;
   if (max<pym)max=pym;
   if (max<pzm)max=pzm;
   if (max<pzm)max=pzm;
   //if tmp is greater or smaller than max and min, set tmp to it.
   //Vp[idx]=RelaxationFactor*(tmp-Vp[idx])+Vp[idx];
   //if need calculate depleted voltage
   double oldP=Vp[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;
   if(tmp<min) {
      Vp[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      Vp[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;
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
void XYZ::GeneralSetup(SquarePointContact &detector)
{
   double dx=detector.Width/(N1-1);
   double dy=detector.Length/(N2-1);
   double dz=detector.Height/(N3-1);
   //general setup
   for(size_t i=0;i<N3;i++)
   {
      for(size_t j=0;j<N2;j++)
      {
         for(size_t k=0;k<N1;k++)
         {
            dC1p.push_back(dx);dC1m.push_back(dx);
            dC2p.push_back(dy);dC2m.push_back(dy);
            dC3p.push_back(dz);dC3m.push_back(dz);
            C1.push_back(k*dx);
            C2.push_back(j*dy);
            C3.push_back(i*dz);
            E1.push_back(0); E2.push_back(0); Et.push_back(0); Vp.push_back(0);
            fIsFixed.push_back(false); fIsDepleted.push_back(false);
            Src.push_back(-detector.GetImpurity(C3[i])*Qe/epsilon);
         }
      }
   }
}
//_____________________________________________________________________________
//
void XYZ::GetInfoFrom(SquarePointContact &spc)
{
   //TODO
   GeneralSetup(spc);
   //boundary
   for(size_t i=0;i<N1*N2*N3;i++)
   {
      if (C1[i]<=0||C1[i]>=spc.Width-1e-5//outer contact
            ||C2[i]<=0||C2[i]>=spc.Length-1e-5
            ||C3[i]>=spc.Height-1e-5)
      {
         fIsDepleted[i]=true;
         fIsFixed[i]=true;
         Vp[i]=spc.Bias[0];
      }
      if(C3[i]>=0&&C3[i]<=spc.PointContactH&&
            C1[i]>=(spc.Width-spc.PointContactW)/2&&
            C1[i]<=(spc.Width+spc.PointContactW)/2&&
            C2[i]>=(spc.Length-spc.PointContactL)/2&&
            C2[i]<=(spc.Length+spc.PointContactL)/2)
      {
         fIsDepleted[i]=true;
         fIsFixed[i]=true;
         Vp[i]=spc.Bias[1];

      }


   }


}
//___________________________________________________________________________
//
void XYZ::CalculateE()
{
   Grid::CalculateE(); // deal with E1
   for (size_t i=0; i<GetN(); i++) { 
      //TODO: E2
      if (i<N1) E2[i]=-(Vp[i+N1]-Vp[i])/dC2p[i]; // lower boundary
      else if (i>GetN()-N1) E2[i]=-(Vp[i]-Vp[i-N1])/dC2m[i]; // upper boundary
      else E2[i]=-(Vp[i+N1]-Vp[i-N1])/(dC2p[i]+dC2m[i]); // the rest E2
      //end E2
      

      //TODO: E1 boundary
      if (i%N1==0) E1[i]=-(Vp[i+1]-Vp[i])/dC1p[i]; // left boundary
      if ((i+1)%N1==0) E1[i]=-(Vp[i]-Vp[i-1])/dC1m[i]; // right boundary
      //end E1
      
      //E3
      if (dC3p[i]==0 || dC3m[i]==0) return;

      if (i<N1*N2) // C3 lower border
         E3[i]=(Vp[i]-Vp[i+N1])/dC3p[i];
      else if (i>=N1*N2*N3-N1*N2) // C3 upper border
         E3[i]=(Vp[i-N1]-Vp[i])/dC3m[i];
      else { // bulk
         E3[i]=(Vp[i-N1]-Vp[i+N1])/(dC3m[i]+dC3p[i]);
      } 
      //end E3
      
      Et[i]=sqrt(E1[i]*E1[i]+E2[i]*E2[i]+E3[i]*E3[i]);//overall E
   }
}
