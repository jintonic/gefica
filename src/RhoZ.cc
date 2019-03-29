#include "RhoZ.h"
#include "Units.h"
#include "PointContact.h"
using namespace GeFiCa;

void RhoZ::GetBoundaryConditionFrom(Detector &detector)
{
   Grid::GetBoundaryConditionFrom(detector); // check number of calls

   TString type(detector.ClassName());
   if (type.Contains("PointContact")==false) {
      Error("GetBoundaryConditionFrom", "%s is not expected. "
            "Please pass in a PointContact detector.", type.Data());
      abort();
   }
   PointContact& pc = (PointContact&) detector;
   pc.CheckConfigurations();
   fDetector = &detector; // for GetC to use fDetector->Bias[]

   for (size_t i=0; i<N1; i++) {
      C1[i]=i*2*pc.Radius/(N1-1);
      dC1p[i]=2*pc.Radius/(N1-1);
      dC1m[i]=2*pc.Radius/(N1-1);
      dC1p.push_back(2*pc.Radius/(N1-1)); dC1m.push_back(2*pc.Radius/(N1-1));
      C1.push_back(pc.BoreR+i*dC1p[i]);
      E1.push_back(0); Et.push_back(0);
      fIsFixed.push_back(false); fIsDepleted.push_back(false);
      Src.push_back(-pc.GetImpurity(C1[i])*Qe/epsilon);
   }
   for (size_t i=0; i<GetN(); i++) {
      if(i>N1-1)C2[i]=C2[i-N1]+pc.Height/(N2-1);
      else C2[i]=0;
      if(i%N1==0)C1[i]=0;
      else C1[i]=C1[i-1]+2*pc.Radius/(N1-1);

      E2[i]=0;
      dC2m[i]=pc.Height/(N2-1);
      dC2p[i]=pc.Height/(N2-1);
   }
   for(size_t i=GetN();i-->0;) C1[i]=C1[i]-pc.Radius;

   dC1m[0]=0; dC1p[N1-1]=0;
   // fix 1st and last points
   fIsFixed[0]=true; fIsFixed[N1-1]=true;
   // linear interpolation between Bias[0] and Bias[1]
   double slope = (pc.Bias[1]-pc.Bias[0])/(N1-1);
   for (size_t i=0; i<N1; i++) Vp.push_back(pc.Bias[0]+slope*i);
   Vp[N1-1]=pc.Bias[1];

   // vacuum
   for (size_t i=0; i<GetN(); i++)
      if (((C1[i]>pc.WrapAroundR-pc.GrooveW && C1[i]<pc.WrapAroundR) ||
               (C1[i]>-pc.WrapAroundR && C1[i]<-pc.WrapAroundR+pc.GrooveW))
            && C2[i]<pc.GrooveH) Src[i]=0;
}
//_____________________________________________________________________________
//
void RhoZ::InitializeGrid()
{
   PointContact& pc = (PointContact&) *fDetector;
   // we want no grid point right on z-axis
   if (N1%2==1) {
      Error("InitializeGrid", "Number of points in D can't be odd! Abort.");
      abort();
   }

   // set initial potential values
   for(size_t i=GetN();i-->0;) {
      Vp[i]=(pc.Bias[1]+pc.Bias[0])/2;
      // set potential for inner electrodes
      if(C1[i]>=-pc.PointContactR && C1[i]<=pc.PointContactR
            && C2[i]<=pc.PointContactH) {
         Vp[i]=pc.Bias[0];
         fIsFixed[i]=true;
      }
   }
   // set potential for outer electrodes
   for(size_t i=GetN()-1;i>=GetN()-N1;i--) {
      fIsFixed[i]=true;
      Vp[i]=pc.Bias[1];
   }
   for(size_t i=0;i<GetN()-N1;i=i+N1) {
      fIsFixed[i]=true;
      fIsFixed[i+N1-1]=true;
      Vp[i]=pc.Bias[1];
      Vp[i+N1-1]=pc.Bias[1];
   }
   for (size_t i=0;i<N1;i++) {
      if(C1[i]>=pc.WrapAroundR||C1[i]<=-pc.WrapAroundR) {
         fIsFixed[i]=true;
         Vp[i]=pc.Bias[1];
      }
   }

   double x1=pc.BoreR+pc.BoreTaperW,
          y1=pc.Height,
          x2=pc.BoreR,
          y2=pc.Height-pc.BoreTaperH,
          x3=pc.Radius-pc.CornerW,
          y3=pc.Height,
          x4=pc.Radius,
          y4=pc.Height-pc.CornerH;
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=(y3-k2*x3);
   if(x3!=x4)
      for(size_t i=0;i<GetN();i++) {
         if(C2[i]>-k2*(C1[i])+b2||(C2[i]>k2*(C1[i])+b2)) {
            fIsFixed[i]=true;
            Vp[i]=pc.Bias[1];
         }
      }
   if(x1!=x2) {
      for (size_t i=0;i<GetN();i++) {
         if(((C2[i]>-k1*(C1[i])+b1 && C2[i]>y2))&&C1[i]<0) {
            fIsFixed[i]=true;
            Vp[i]=pc.Bias[1];
         }
         if(((C2[i]>k1*(C1[i])+b1 && C2[i]>y2))&&C1[i]>0) {
            fIsFixed[i]=true;
            Vp[i]=pc.Bias[1];
         }
      }
   } else { //x1==x2
      for (size_t i=0;i<GetN();i++) {
         if(C1[i]<=x1&&C1[i]>=-x1&&C2[i]>=y2) {
            fIsFixed[i]=true;
            Vp[i]=pc.Bias[1];
         }
      }
   }
   for (size_t i=0;i<GetN();i++) {
      if(C1[i]<=pc.BoreR&&C1[i]>=-pc.BoreR&&C2[i]>=pc.Height-pc.BoreH) {
         fIsFixed[i]=true;
         Vp[i]=pc.Bias[1];
      }
   }
   SetBoundary();
}
//_____________________________________________________________________________
//
void RhoZ::OverRelaxAt(size_t idx)
{
   PointContact& pc = (PointContact&) *fDetector;
   if (fIsFixed[idx])return; 
   // 2nd-order Successive Over-Relaxation
   double density=Src[idx];
   double drm=dC1m[idx]; // dr_minus
   double drp=dC1p[idx];
   double dzm=dC2m[idx];
   double dzp=dC2p[idx];
   double pzm,pzp,prm,prp; // pzm: potential_z_plus
   if(idx>=N1)pzm=Vp[idx-N1];
   else pzm=Vp[idx+N1];
   if(idx>=GetN()-N1)pzp=Vp[idx];
   else pzp=Vp[idx+N1];
   if(idx%N1==0)prm=Vp[idx];
   else prm=Vp[idx-1];
   if(idx%N1==N1-1)prp=Vp[idx];
   else prp=Vp[idx+1];
   double tmp=(density/epsilon
         + 1/C1[idx]*(prp-prm)/(drm+drp) +(prp/drp+prm/drm)*2/(drm+drp)
         + (pzp/dzp+pzm/dzm)*2/(dzp+dzm))/
      ((1/drm+1/drp)*2/(drm+drp)+(1/dzp+1/dzm)*2/(dzp+dzm));
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=prm;
   double max=prm;
   if(min>prp)min=prp;
   if (min>pzp)min=pzp;
   if (min>pzm)min=pzm;

   //find max
   if(max<prp)max=prp;
   if (max<pzp)max=pzp;
   if (max<pzm)max=pzm;
   //if tmp is greater or smaller than max and min, set tmp to it.

   //over relax
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

   if(fIsDepleted[idx]||pc.Bias[0]==pc.Bias[1]) Vp[idx]=tmp;
}
//_____________________________________________________________________________
//
double RhoZ::GetC()
{
   Grid::GetC(); // calculate field excluding undepleted region

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   PointContact& pc = (PointContact&) *fDetector;
   double dV=pc.Bias[0]-pc.Bias[1]; if(dV<0)dV=-dV;
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
//______________________________________________________________________________
//
void RhoZ::SetBoundary()
{
   PointContact& pc = (PointContact&) *fDetector;
   for(size_t i=0;i<GetN();i++) {
      if (C2[i]-pc.PointContactH<dC2m[i]&&C2[i]>pc.PointContactH
            &&C1[i]<pc.PointContactR&&C1[i]>-pc.PointContactR)
         dC2m[i]=C2[i]-pc.PointContactH;
      if(C1[i]-pc.PointContactR<dC1m[i]&&C1[i]>0&&C2[i]<pc.PointContactH)
         dC1m[i]=C1[i]-pc.PointContactR;
      if(-C1[i]-pc.PointContactR<dC1p[i]&&C1[i]<0&&C2[i]<pc.PointContactH)
         dC1p[i]=-C1[i]-pc.PointContactR;
      if(pc.WrapAroundR-C1[i]<dC1p[i]&&C1[i]<pc.WrapAroundR&&i<N1)
         dC1p[i]=pc.WrapAroundR-C1[i];
      if(pc.WrapAroundR+C1[i]<dC1p[i]&&C1[i]>-pc.WrapAroundR&&i<N1)
         dC1m[i]=pc.WrapAroundR+C1[i];
   }
   double k=pc.TaperH/(pc.TaperW);
   double b=-(pc.Radius-pc.TaperW)*k;

   for(size_t i=0;i<GetN();i++) {
      if(C2[i]<=C1[i]*k+b) {
         fIsFixed[i]=true;
         Vp[i]=pc.Bias[1];
      }
      if(C2[i]<=-C1[i]*k+b) {
         fIsFixed[i]=true;
         Vp[i]=pc.Bias[1];
      }
      if(C2[i]-(C1[i]*k+b)<dC2p[i]) {
         dC2m[i]=C2[i]-(k*C1[i]+b);
         dC1p[i]=C1[i]-b/k-C2[i]/k;
      }
      if(C2[i]-(-k*C1[i]+b)<dC2m[i]) {
         dC2m[i]=C2[i]-(-C1[i]*k+b);
         dC1m[i]=-C1[i]/k-b/k-C2[i];
      }
   }
   double x1=pc.BoreTaperW,
          y1=pc.Height,
          x2=pc.BoreR,
          y2=pc.Height-pc.BoreTaperH,
          x3=pc.Radius-pc.CornerW,
          y3=pc.Height,
          x4=pc.Radius,
          y4=pc.Height-pc.CornerH;
   // y = k x + b
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=y3-k2*x3;

   for (size_t i=0;i<GetN();i++) {
      if (x1!=x2) {
         //right side of hole taper
         if (C1[i]-C2[i]/k1+b1/k1<dC1m[i] && C2[i]>y2 &&
               C1[i]-C2[i]/k1+b1/k1>0)
            dC1m[i]=C1[i]-C2[i]/k1+b1/k1;
         //left side of hole taper
         if (-C1[i]-C2[i]/k1+b1/k1>0
               &&-C1[i]-C2[i]/k1+b1/k1<dC1p[i]&&C2[i]>y2)
            dC1p[i]=-C1[i]-C2[i]/k1+b1/k1;
      } else { //x1==x2
         //right side of hole taper
         if (C1[i]-x1<dC1m[i] && C2[i]>y2 && C1[i]-x1>0) dC1m[i]=C1[i]-x1;
         //left side of hole taper
         if (-C1[i]-x1>0&&-C1[i]-x1<dC1p[i]&&C2[i]>y2) dC1p[i]=-C1[i]-x1;
      }
      //right side of hole taper
      if (C1[i]-pc.BoreR<dC1m[i] && C2[i]>pc.BoreH && C1[i]-pc.BoreR>0)
         dC1m[i]=C1[i]-pc.BoreR;
      //left side of hole
      if (-C1[i]-pc.BoreR>0&&-C1[i]-pc.BoreR<dC1p[i]&&C2[i]>pc.BoreH)
         dC1p[i]=-C1[i]-pc.BoreR;
      //left corner
      if (C1[i]+C2[i]/k2-b2/k2>0 && C1[i]+C2[i]/k2-b2/k2<dC1m[i] &&
            C2[i]>y4) dC1m[i]=C1[i]+C2[i]/k2-b2/k2;
      //right corner
      if (-C1[i]+C2[i]/k2-b2/k2>0
            &&-C1[i]+C2[i]/k2-b2/k2<dC1p[i]&&C2[i]>y4)
         dC1p[i]=-C1[i]+C2[i]/k2-b2/k2;
      //down right side of hole
      if (-C2[i]+C1[i]*k1+b1>0&&-C2[i]+C1[i]*k1+b1<dC2p[i]&&C2[i]>y2)
         dC2p[i]=-C2[i]+C1[i]*k1+b1;
      //down right of corner
      if (-C2[i]-C1[i]*k2+b2>0&&-C2[i]-C1[i]*k2+b2<dC2p[i]&&C2[i]>y4)
         dC2p[i]=-C2[i]-C1[i]*k2+b2;
      //down left side of hole
      if (-C2[i]-C1[i]*k1+b1>0&&-C2[i]-C1[i]*k1+b1<dC2p[i]&&C2[i]>y2)
         dC2p[i]=-C2[i]-C1[i]*k1+b1;
      //down left of corner
      if (-C2[i]+C1[i]*k2+b2>0&&-C2[i]+C1[i]*k2+b2<dC2p[i]&&C2[i]>y4)
         dC2p[i]=-C2[i]+C1[i]*k2+b2;
      //down center of hole
      if (y2-C2[i]<dC2p[i]&&C1[i]>-pc.BoreR&&C1[i]<pc.BoreR)
         dC2p[i]=y2-C2[i];
   }
}
//_____________________________________________________________________________
//

