#include <stdlib.h>

#include "PointContactDZ.h"
#include "iostream"
#include "Units.h"
using namespace GeFiCa;
using namespace std;
//for the case which point contact may not fall right on the grid, it will case
//a great different potential on space close to point contact the current
//design does not include when point contact is close to boundary of the whole
//detector.
PointContactDZ::PointContactDZ(int nd, int nz, const char *name,
      const char *title) : RhoZ(nd, nz, name, title),
   Height(5*cm),
   Radius(3*cm),
   PointContactH(0.01*cm),
   PointContactR(0.1*cm),
   HoleH(0),
   HoleInnerR(0),
   HoleOuterR(0),
   TaperW(0.3*cm),
   TaperH(0.3*cm),
   CornerW(0.3*cm),
   CornerH(0.3*cm),
   WrapArroundR(0.5*cm) {};
//_____________________________________________________________________________
//
void PointContactDZ::SetBoundary()
{
   for(int i=0;i<n;i++)
   {
      if(fC2[i]-PointContactH<fdC2m[i]&&fC2[i]>PointContactH&&fC1[i]<PointContactR&&fC1[i]>-PointContactR)
      {
         fdC2m[i]=fC2[i]-PointContactH;
      }
      if(fC1[i]-PointContactR<fdC1m[i]&&fC1[i]>0&&fC2[i]<PointContactH)
      {
         fdC1m[i]=fC1[i]-PointContactR;
      }
      if(-fC1[i]-PointContactR<fdC1p[i]&&fC1[i]<0&&fC2[i]<PointContactH)
      {
         fdC1p[i]=-fC1[i]-PointContactR;
      }
      if(WrapArroundR-fC1[i]<fdC1p[i]&&fC1[i]<WrapArroundR)
      {
         fdC1p[i]=WrapArroundR-fC1[i];
      }
      if(WrapArroundR+fC1[i]<fdC1p[i]&&fC1[i]>-WrapArroundR)
      {
         fdC1p[i]=WrapArroundR+fC1[i];
      }
   }
   double k=TaperH/(TaperW);
   double b=-(Radius-TaperW)*k;

   for(int i=0;i<n;i++)
   {
      if(fC2[i]<=fC1[i]*k+b)
      {
         fIsFixed[i]=true;
         fV[i]=V0;
      }
      if(fC2[i]<=-fC1[i]*k+b)
      {
         fIsFixed[i]=true;
         fV[i]=V0;
      }
      if(fC2[i]-(fC1[i]*k+b)<fdC2m[i])
      {
         fdC2m[i]=fC2[i]-(k*fC1[i]+b);
         fdC1p[i]=fC2[i]-b-k*fC1[i];
      }
      if(fC2[i]-(-k*fC1[i]+b)<fdC2m[i])
      {
         fdC2m[i]=fC2[i]-(-fC1[i]*k+b);
         fdC1m[i]=-fC2[i]*k+b-fC1[i];
      }
   }
   double x1=HoleOuterR,
          y1=Height,
          x2=HoleInnerR,
          y2=Height-HoleH,
          x3=Radius-CornerW,
          y3=Height,
          x4=Radius,
          y4=Height-CornerH;
   // y = k x + b
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=y3-k2*x3;

   for (int i=0;i<n;i++) {
      //right side of hole
      if(fC1[i]-fC2[i]/k1+b1/k1<fdC1m[i] && fC2[i]>y2 &&
            fC1[i]-fC2[i]/k1+b1/k1>0) fdC1m[i]=fC1[i]-fC2[i]/k1+b1/k1;
      //left corner
      if(fC1[i]+fC2[i]/k2-b2/k2>0 && fC1[i]+fC2[i]/k2-b2/k2<fdC1m[i] &&
            fC2[i]>y4) fdC1m[i]=fC1[i]+fC2[i]/k2-b2/k2;
      //left side of hole
      if(-fC1[i]-fC2[i]/k1+b1/k1>0&&-fC1[i]-fC2[i]/k1+b1/k1<fdC1p[i]&&fC2[i]>y2)
         fdC1p[i]=-fC1[i]-fC2[i]/k1+b1/k1;
      //right corner
      if(-fC1[i]+fC2[i]/k2-b2/k2>0&&-fC1[i]+fC2[i]/k2-b2/k2<fdC1p[i]&&fC2[i]>y4)
         fdC1p[i]=-fC1[i]+fC2[i]/k2-b2/k2;
      //down right side of hole
      if(-fC2[i]+fC1[i]*k1+b1>0&&-fC2[i]+fC1[i]*k1+b1<fdC2p[i]&&fC2[i]>y2)
         fdC2p[i]=-fC2[i]+fC1[i]*k1+b1;
      //down right of corner
      if(-fC2[i]-fC1[i]*k2+b2>0&&-fC2[i]-fC1[i]*k2+b2<fdC2p[i]&&fC2[i]>y4)
         fdC2p[i]=-fC2[i]-fC1[i]*k2+b2;
      //down left side of hole
      if(-fC2[i]-fC1[i]*k1+b1>0&&-fC2[i]-fC1[i]*k1+b1<fdC2p[i]&&fC2[i]>y2)
         fdC2p[i]=-fC2[i]-fC1[i]*k1+b1;
      //down left of corner
      if(-fC2[i]+fC1[i]*k2+b2>0&&-fC2[i]+fC1[i]*k2+b2<fdC2p[i]&&fC2[i]>y4)
         fdC2p[i]=-fC2[i]+fC1[i]*k2+b2;
      //down center of hole
      if(y2-fC2[i]<fdC2p[i]&&fC1[i]>-HoleInnerR&&fC1[i]<HoleInnerR)
         fdC2p[i]=y2-fC2[i];
   }
}
//_____________________________________________________________________________
//
void PointContactDZ::Initialize()
{
   // we want no grid point right on z-axis
   if (n1%2==1) Fatal("Initialize", "Number of points in D cannot be odd!");

   SetStepLength(2*Radius/(n1-1),Height/(n2-1));
   for(int i=n;i-->0;) fC1[i]=fC1[i]-Radius;

   // set initial potential values
   for(int i=n;i-->0;) {
      fV[i]=(V0+V1)/2;
      // set potential for inner electrodes
      if(fC1[i]>=-PointContactR && fC1[i]<=RointContactR
            && fC2[i]<=PointContactH) {
         fV[i]=V1;
         fIsFixed[i]=true;
      }
   }
   // set potential for outer electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fV[i]=V0;
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fV[i]=V0;
      fV[i+n1-1]=V0;
   }
   for (int i=0;i<n1;i++) {
      if(fC1[i]>=WrapArroundR||fC1[i]<=-WrapArroundR) {
         fIsFixed[i]=true;
         fV[i]=V0;
      }
   }

   double x1=HoleOuterR,
          y1=Height,
          x2=HoleInnerR,
          y2=Height-HoleH,
          x3=Radius-CornerW,
          y3=Height,
          x4=Radius,
          y4=Height-CornerH;
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=(y3-k2*x3);

   for (int i=0;i<n;i++) {
      if(((fC2[i]>-k1*(fC1[i])+b1 
                  && fC2[i]>y2)||(fC2[i]>-k2*(fC1[i])+b2))&&fC1[i]<0) {
         fIsFixed[i]=true;
         fV[i]=V0;
      }
      if(((fC2[i]>k1*(fC1[i])+b1
                  && fC2[i]>y2)||(fC2[i]>k2*(fC1[i])+b2))&&fC1[i]>0) {
         fIsFixed[i]=true;
         fV[i]=V0;
      }

   }
   SetBoundary();
}
//_____________________________________________________________________________
//
bool PointContactDZ::CalculateField(int idx)
{
   if (!XY::CalculateField(idx)) return false;

   if (fC2[idx]>PointContactH-fdC2m[idx]
         && fC2[idx]<PointContactH+fdC2p[idx]) // PC top boundary
      fE2[idx]=(fV[idx]-fV[idx+n1])/fdC2p[idx];
   if (fC1[idx]>-PointContactR-fdC1m[idx]
         && fC1[idx]<-PointContactR+fdC1p[idx]) // PC left boundary
      fE1[idx]=(fV[idx]-fV[idx-1])/fdC1m[idx];
   if (fC1[idx]>PointContactR-fdC1m[idx]
         && fC1[idx]<PointContactR+fdC1p[idx]) // PC right boundary
      fE1[idx]=(fV[idx]-fV[idx+1])/fdC1p[idx];

   return true;
}
//_____________________________________________________________________________
//
#include <fstream>
bool PointContactDZ::SaveFieldAsFieldgen(const char * fout)
{
   ofstream outfile(fout);

   outfile<<"# height "<< Height;        
   outfile<<"\n# xtal_radius "<<Radius;
   outfile<<"\n# pc_length   "<<PointContactH;        
   outfile<<"\n# pc_radius   "<<PointContactR;         
   outfile<<"\n# wrap_around_radius "<<WrapArroundR; 
   outfile<<"\n# grid size on r "<<fdC1p[0];
   outfile<<"\n# grid size on z "<<fdC2p[0];
   outfile<<"\n# impurity_z0  "<<fImpurity[0];
   outfile<<"\n# xtal_HV      "<<V1;
   outfile<<"\n# max_iterations "<<MaxIterations;
   outfile<<"\n# ";
   outfile<<"\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)";
   for (int i=0;i<n;i++) {
      double E=sqrt(fE1[i]*fE1[i]+fE2[i]*fE2[i]);
      outfile<<"\n"<<fC1[i]<<"  "<<fC2[i]<<"  "<<fV[i]<<"  "<<E<<"  "<<fE1[i]<<"  "<<fE2[i];
   }
   outfile.close();
   return true;
}
