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
   HoleR(0), 
   HoleTaperW(0),
   HoleTaperH(0),
   TaperW(0.3*cm),
   TaperH(0.3*cm),
   CornerW(0.3*cm),
   CornerH(0.5*cm),
   WrapArroundR(2.5*cm),
   GrooveW(0.1*cm), 
   GrooveH(0.1){};
//_____________________________________________________________________________
//
void PointContactDZ::SetBoundary()
{
   for(int i=0;i<fN;i++)
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
      if(WrapArroundR-fC1[i]<fdC1p[i]&&fC1[i]<WrapArroundR&&i<fN1)
      {
         fdC1p[i]=WrapArroundR-fC1[i];
      }
      if(WrapArroundR+fC1[i]<fdC1p[i]&&fC1[i]>-WrapArroundR&&i<fN1)
      {
         fdC1m[i]=WrapArroundR+fC1[i];
      }
   }
   double k=TaperH/(TaperW);
   double b=-(Radius-TaperW)*k;

   for(int i=0;i<fN;i++)
   {
      if(fC2[i]<=fC1[i]*k+b)
      {
         fIsFixed[i]=true;
         fV[i]=V1;
      }
      if(fC2[i]<=-fC1[i]*k+b)
      {
         fIsFixed[i]=true;
         fV[i]=V1;
      }
      if(fC2[i]-(fC1[i]*k+b)<fdC2p[i])
      {
         fdC2m[i]=fC2[i]-(k*fC1[i]+b);
         fdC1p[i]=fC1[i]-b/k-fC2[i]/k;
      }
      if(fC2[i]-(-k*fC1[i]+b)<fdC2m[i])
      {
         fdC2m[i]=fC2[i]-(-fC1[i]*k+b);
         fdC1m[i]=-fC1[i]/k-b/k-fC2[i];
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
      if (x1!=x2)
      {
         //right side of hole taper
         if(fC1[i]-fC2[i]/k1+b1/k1<fdC1m[i] && fC2[i]>y2 &&
               fC1[i]-fC2[i]/k1+b1/k1>0)
            fdC1m[i]=fC1[i]-fC2[i]/k1+b1/k1;
         //left side of hole taper
         if(-fC1[i]-fC2[i]/k1+b1/k1>0&&-fC1[i]-fC2[i]/k1+b1/k1<fdC1p[i]&&fC2[i]>y2)
            fdC1p[i]=-fC1[i]-fC2[i]/k1+b1/k1;
      }
      else//x1==x2
      {
         //right side of hole taper
         if(fC1[i]-x1<fdC1m[i] && fC2[i]>y2 &&
               fC1[i]-x1>0)
            fdC1m[i]=fC1[i]-x1;
         //left side of hole taper
         if(-fC1[i]-x1>0&&-fC1[i]-x1<fdC1p[i]&&fC2[i]>y2)
            fdC1p[i]=-fC1[i]-x1;
      }
      //right side of hole taper
      if(fC1[i]-HoleR<fdC1m[i] && fC2[i]>HoleH &&
            fC1[i]-HoleR>0)
          fdC1m[i]=fC1[i]-HoleR;
      //left side of hole
      if(-fC1[i]-HoleR>0&&-fC1[i]-HoleR<fdC1p[i]&&fC2[i]>HoleH)
         fdC1p[i]=-fC1[i]-HoleR;
      //left corner
      if(fC1[i]+fC2[i]/k2-b2/k2>0 && fC1[i]+fC2[i]/k2-b2/k2<fdC1m[i] &&
            fC2[i]>y4) fdC1m[i]=fC1[i]+fC2[i]/k2-b2/k2;
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
      if(y2-fC2[i]<fdC2p[i]&&fC1[i]>-HoleR&&fC1[i]<HoleR)
         fdC2p[i]=y2-fC2[i];
   }
}
//_____________________________________________________________________________
//
void PointContactDZ::InitializeGrid()
{
   // we want no grid point right on z-axis
   if (fN1%2==1) Fatal("InitializeGrid", "Number of points in D can't be odd!");

   SetStepLength(2*Radius/(fN1-1),Height/(fN2-1));
   for(int i=fN;i-->0;) fC1[i]=fC1[i]-Radius;

   // set initial potential values
   for(int i=fN;i-->0;) {
      fV[i]=(V1+V0)/2;
      // set potential for inner electrodes
      if(fC1[i]>=-PointContactR && fC1[i]<=PointContactR
            && fC2[i]<=PointContactH) {
         fV[i]=V0;
         fIsFixed[i]=true;
      }
   }
   // set potential for outer electrodes
   for(int i=fN-1;i>=fN-fN1;i--) {
      fIsFixed[i]=true;
      fV[i]=V1;
   }
   for(int i=0;i<fN-fN1;i=i+fN1) {
      fIsFixed[i]=true;
      fIsFixed[i+fN1-1]=true;
      fV[i]=V1;
      fV[i+fN1-1]=V1;
   }
   for (int i=0;i<fN1;i++) {
      if(fC1[i]>=WrapArroundR||fC1[i]<=-WrapArroundR) {
         fIsFixed[i]=true;
         fV[i]=V1;
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
      for(int i=0;i<fN;i++)
      {
         if(fC2[i]>-k2*(fC1[i])+b2||(fC2[i]>k2*(fC1[i])+b2))
         {
            fIsFixed[i]=true;
            fV[i]=V1;
         }
      }
   if(x1!=x2)
   {
      for (int i=0;i<fN;i++) {
         if(((fC2[i]>-k1*(fC1[i])+b1 
                     && fC2[i]>y2))&&fC1[i]<0) {
            fIsFixed[i]=true;
            fV[i]=V1;
         }
         if(((fC2[i]>k1*(fC1[i])+b1
                     && fC2[i]>y2))&&fC1[i]>0) {
            fIsFixed[i]=true;
            fV[i]=V1;
         }
      }
   }
   else//x1==x2
   {
      for (int i=0;i<fN;i++) 
      {
         if(fC1[i]<=x1&&fC1[i]>=-x1&&fC2[i]>=y2)
         {
            fIsFixed[i]=true;
            fV[i]=V1;
         }
      }

   }
   for (int i=0;i<fN;i++) 
   {
      if(fC1[i]<=HoleR&&fC1[i]>=-HoleR&&fC2[i]>=Height-HoleH)
      {
         fIsFixed[i]=true;
         fV[i]=V1;
      }
   }
   SetBoundary();
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
   outfile<<"\n# xtal_HV      "<<V0;
   outfile<<"\n# max_iterations "<<MaxIterations;
   outfile<<"\n# ";
   outfile<<"\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)";
   for (int i=0;i<fN;i++) {
      double E=sqrt(fE1[i]*fE1[i]+fE2[i]*fE2[i]);
      outfile<<"\n"<<fC1[i]<<"  "<<fC2[i]<<"  "<<fV[i]<<"  "<<E<<"  "<<fE1[i]<<"  "<<fE2[i];
   }
   outfile.close();
   return true;
}
void PointContactDZ::SetGridImpurity()
{

   X::SetGridImpurity();
   
   if(TaperW+WrapArroundR<Radius)
   {
      for(int i=0;i<fN;i++)
      {
         if((fC1[i]>WrapArroundR-GrooveW && fC1[i]<WrapArroundR && fC2[i]<GrooveH)||
            (fC1[i]<-(WrapArroundR-GrooveW) && fC1[i]>-(WrapArroundR) && fC2[i]<GrooveH))
         {
            fImpurity[i]=0;
            fV[i]=0;
         }
      }
   }
}
