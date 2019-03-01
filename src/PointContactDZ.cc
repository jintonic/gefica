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
//_____________________________________________________________________________
//
void PointContactDZ::SetBoundary()
{
   for(int i=0;i<n;i++)
   {
      if(fC2[i]-PointContactZ<fdC2m[i]&&fC2[i]>PointContactZ&&fC1[i]<PointContactR&&fC1[i]>-PointContactR)
      {
         fdC2m[i]=fC2[i]-PointContactZ;
      }
      if(fC1[i]-PointContactR<fdC1m[i]&&fC1[i]>0&&fC2[i]<PointContactZ)
      {
         fdC1m[i]=fC1[i]-PointContactR;
      }
      if(-fC1[i]-PointContactR<fdC1p[i]&&fC1[i]<0&&fC2[i]<PointContactZ)
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
   double k=TaperZ/(TaperLength);
   double b=-(Radius-TaperLength)*k;

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
   double x1=HoleOutterR,
	  y1=Z,
	  x2=HoleInnerR,
	  y2=Z-HoleZ,
	  x3=Radius-ConnorLength,
	  y3=Z,
	  x4=Radius,
	  y4=Z-ConnorZ;
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
   if (n1%2==1) Fatal("Initialize", "Number of grids in D cannot be odd!");

   if (Z0>=Z) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            Z0, Z);
      return;
   }
   double RUpperBound,RLowerBound,PointBegin,PointEnd;
   RUpperBound=Radius;
   RLowerBound=-Radius;
   PointBegin=-PointContactR;
   PointEnd=PointContactR;
   double steplength1=(RUpperBound-RLowerBound)/(n1-1);
   double steplength2=(Z-Z0)/(n2-1);
   SetStepLength(steplength1,steplength2);
   for(int i=n;i-->0;) 
   {
      fC1[i]=fC1[i]+RLowerBound;
   } 

   // set initial potential values
   for(int i=n;i-->0;) {
      fV[i]=(V0+V1)/2;//common this line for finding depleat voltage
      // set potential for inner electrodes
      if(fC1[i]>=PointBegin&&fC1[i]<=PointEnd&&fC2[i]<=PointContactZ) {
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
   for (int i=0;i<n1;i++)
   {
      if(fC1[i]>=WrapArroundR||fC1[i]<=-WrapArroundR)
      {
         fIsFixed[i]=true;
         fV[i]=V0;
      }
   }

   double x1=HoleOutterR,
	  y1=Z,
	  x2=HoleInnerR,
	  y2=Z-HoleZ,
	  x3=Radius-ConnorLength,
	  y3=Z,
	  x4=Radius,
	  y4=Z-ConnorZ;
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=(y3-k2*x3);

   for (int i=0;i<n;i++)
   {
     if(((fC2[i]>-k1*(fC1[i])+b1  && fC2[i]>y2)||(fC2[i]>-k2*(fC1[i])+b2))&&fC1[i]<0)
     {
       fIsFixed[i]=true;
       fV[i]=V0;
     }
     if(((fC2[i]>k1*(fC1[i])+b1  && fC2[i]>y2)||(fC2[i]>k2*(fC1[i])+b2))&&fC1[i]>0)
     {
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

   if (fC2[idx]>PointContactZ-fdC2m[idx]
         && fC2[idx]<PointContactZ+fdC2p[idx]) // PC top boundary
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

   outfile<<"# height "<< Z-Z0;        
   outfile<<"\n# xtal_radius "<<Radius;
   outfile<<"\n# pc_length   "<<PointContactZ;        
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
