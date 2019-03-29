#include "Units.h"
#include "PointContact.h"
using namespace GeFiCa;
//for the case which point contact may not fall right on the grid, it will case
//a great different potential on space close to point contact the current
//design does not include when point contact is close to boundary of the whole
//detector.
PointContactDZ::PointContactDZ(int nd, int nz, const char *name,
      const char *title) : RhoZ(nd, nz, name, title),
   Height(5*cm),
   Radius(3*cm),
   PointContactH(0),
   PointContactR(1*mm),
   HoleH(0),
   HoleR(0), 
   HoleTaperW(0),
   HoleTaperH(0),
   TaperW(1*mm),
   TaperH(1*mm),
   CornerW(1*mm),
   CornerH(1*mm),
   WrapAroundR(-1),
   GrooveW(0), 
   GrooveH(0) {}
//_____________________________________________________________________________
//
void PointContactDZ::SetBoundary()
{
   for(int i=0;i<fN;i++) {
      if (C2[i]-PointContactH<dC2m[i]&&C2[i]>PointContactH
            &&C1[i]<PointContactR&&C1[i]>-PointContactR)
         dC2m[i]=C2[i]-PointContactH;
      if(C1[i]-PointContactR<dC1m[i]&&C1[i]>0&&C2[i]<PointContactH)
         dC1m[i]=C1[i]-PointContactR;
      if(-C1[i]-PointContactR<dC1p[i]&&C1[i]<0&&C2[i]<PointContactH)
         dC1p[i]=-C1[i]-PointContactR;
      if(WrapAroundR-C1[i]<dC1p[i]&&C1[i]<WrapAroundR&&i<N1)
         dC1p[i]=WrapAroundR-C1[i];
      if(WrapAroundR+C1[i]<dC1p[i]&&C1[i]>-WrapAroundR&&i<N1)
         dC1m[i]=WrapAroundR+C1[i];
   }
   double k=TaperH/(TaperW);
   double b=-(Radius-TaperW)*k;

   for(int i=0;i<fN;i++) {
      if(C2[i]<=C1[i]*k+b) {
         fIsFixed[i]=true;
         V[i]=Bias[1];
      }
      if(C2[i]<=-C1[i]*k+b) {
         fIsFixed[i]=true;
         V[i]=Bias[1];
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
   double x1=HoleTaperW,
          y1=Height,
          x2=HoleR,
          y2=Height-HoleTaperH,
          x3=Radius-CornerW,
          y3=Height,
          x4=Radius,
          y4=Height-CornerH;
   // y = k x + b
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=y3-k2*x3;

   for (int i=0;i<fN;i++) {
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
      if (C1[i]-HoleR<dC1m[i] && C2[i]>HoleH && C1[i]-HoleR>0)
          dC1m[i]=C1[i]-HoleR;
      //left side of hole
      if (-C1[i]-HoleR>0&&-C1[i]-HoleR<dC1p[i]&&C2[i]>HoleH)
         dC1p[i]=-C1[i]-HoleR;
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
      if (y2-C2[i]<dC2p[i]&&C1[i]>-HoleR&&C1[i]<HoleR)
         dC2p[i]=y2-C2[i];
   }
}
//_____________________________________________________________________________
//
void PointContactDZ::InitializeGrid()
{
   // we want no grid point right on z-axis
   if (N1%2==1) {
      Error("InitializeGrid", "Number of points in D can't be odd! Abort.");
      abort();
   }

   if (WrapAroundR<0) WrapAroundR=Radius-TaperW;
   SetStepLength(2*Radius/(N1-1),Height/(N2-1));
   for(int i=fN;i-->0;) C1[i]=C1[i]-Radius;

   // set initial potential values
   for(int i=fN;i-->0;) {
      V[i]=(Bias[1]+Bias[0])/2;
      // set potential for inner electrodes
      if(C1[i]>=-PointContactR && C1[i]<=PointContactR
            && C2[i]<=PointContactH) {
         V[i]=Bias[0];
         fIsFixed[i]=true;
      }
   }
   // set potential for outer electrodes
   for(int i=fN-1;i>=fN-N1;i--) {
      fIsFixed[i]=true;
      V[i]=Bias[1];
   }
   for(int i=0;i<fN-N1;i=i+N1) {
      fIsFixed[i]=true;
      fIsFixed[i+N1-1]=true;
      V[i]=Bias[1];
      V[i+N1-1]=Bias[1];
   }
   for (int i=0;i<N1;i++) {
      if(C1[i]>=WrapAroundR||C1[i]<=-WrapAroundR) {
         fIsFixed[i]=true;
         V[i]=Bias[1];
      }
   }

   double x1=HoleR+HoleTaperW,
          y1=Height,
          x2=HoleR,
          y2=Height-HoleTaperH,
          x3=Radius-CornerW,
          y3=Height,
          x4=Radius,
          y4=Height-CornerH;
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=(y3-k2*x3);
   if(x3!=x4)
      for(int i=0;i<fN;i++) {
         if(C2[i]>-k2*(C1[i])+b2||(C2[i]>k2*(C1[i])+b2)) {
            fIsFixed[i]=true;
            V[i]=Bias[1];
         }
      }
   if(x1!=x2) {
      for (int i=0;i<fN;i++) {
         if(((C2[i]>-k1*(C1[i])+b1 && C2[i]>y2))&&C1[i]<0) {
            fIsFixed[i]=true;
            V[i]=Bias[1];
         }
         if(((C2[i]>k1*(C1[i])+b1 && C2[i]>y2))&&C1[i]>0) {
            fIsFixed[i]=true;
            V[i]=Bias[1];
         }
      }
   } else { //x1==x2
      for (int i=0;i<fN;i++) {
         if(C1[i]<=x1&&C1[i]>=-x1&&C2[i]>=y2) {
            fIsFixed[i]=true;
            V[i]=Bias[1];
         }
      }
   }
   for (int i=0;i<fN;i++) {
      if(C1[i]<=HoleR&&C1[i]>=-HoleR&&C2[i]>=Height-HoleH) {
         fIsFixed[i]=true;
         V[i]=Bias[1];
      }
   }
   SetBoundary();
}
//_____________________________________________________________________________
//
#include <cmath>
#include <fstream>
using namespace std;
void PointContactDZ::Export2fieldgen(const char *output)
{
   ofstream file(output);
   file<<"# xtal_length "<<Height/mm<<endl;
   file<<"# xtal_radius "<<Radius/mm<<endl;
   file<<"# pc_length   "<<PointContactH/mm<<endl;
   file<<"# pc_radius   "<<PointContactR/mm<<endl;
   file<<"# taper_length "<<TaperH/mm<<endl;
   file<<"# wrap_around_radius "<<WrapAroundR/mm<<endl;
   file<<"# ditch_depth "<<GrooveH/mm<<endl;
   file<<"# ditch_thickness "<<GrooveW/mm<<endl;
   file<<"# xtal_grid "<<dC1p[0]/mm<<endl;
   file<<"# impurity_z0 "<<fImpurity[0]*cm3<<endl;
   file<<"# impurity_gradient "
      <<(fImpurity[fN-1]-fImpurity[0])/dC2p[0]*cm*cm3<<endl;
   file<<"# xtal_HV      "<<abs(Bias[0]-Bias[1])/volt<<endl;
   file<<"# max_iterations "<<GetNsor()<<endl;
   file<<"#"<<endl;
   file<<"## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)"<<endl;
   for (int i=0;i<fN;i++) file<<C1[i]<<" "<<C2[i]<<" "<<V[i]<<" "
      <<sqrt(E1[i]*E1[i]+E2[i]*E2[i])<<" "<<E1[i]<<" "<<E2[i]<<endl;
   file.close();
}
//_____________________________________________________________________________
//
void PointContactDZ::SetGridImpurity()
{
   X::SetGridImpurity();

   if (TaperW+WrapAroundR>Radius) {
      Error("SetGridImpurity", "TaperW(%.1fmm) + WrapAroundR(%.1fmm)"
            " > Radiu(%.1fmm). Abort!", TaperW/mm, WrapAroundR/mm, Radius/mm);
      abort();
   } else if (WrapAroundR-GrooveW<PointContactR) {
      Error("SetGridImpurity", "WrapAroundR(%.1fmm) - GrooveW(%.1fmm)"
            " < PointContactR(%.1fmm). Abort!",
            WrapAroundR/mm, GrooveW/mm, PointContactR/mm);
      abort();
   } else {
      for (int i=0; i<fN; i++) {
         if (((C1[i]>WrapAroundR-GrooveW && C1[i]<WrapAroundR) ||
                  (C1[i]>-WrapAroundR && C1[i]<-WrapAroundR+GrooveW))
               && C2[i]<GrooveH) fImpurity[i]=0;
      }
   }
}
//_____________________________________________________________________________
//
double PointContactDZ::GetData(double x, double y, double z, double *data)
{
   //if (point in boundary) return Bias[1];
   return XY::GetData(x,y,z,data);
}
