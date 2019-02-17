#include "ReversedCoaxialRZ.h"
#include "iostream"
#include "Units.h"
#include <cmath>
using namespace GeFiCa;
void ReversedCoaxialRZ::SetupBoundary()
{
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
void ReversedCoaxialRZ::Initialize()
{
   if (Radius<=HoleOutterR||Radius<=HoleInnerR) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            Radius, HoleOutterR);
      return;
   }
   double steplength1=(Radius*2)/(n1-1);
   double steplength2=(Z-Z0)/(n2-1);
   SetStepLength(steplength1,steplength2);
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

   for(int i=n;i-->0;) 
   {
      fC1[i]=fC1[i]-Radius;
      fV[i]=(V0+V1)/2;
   }
   // set potential for electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fV[i]=V0;
      if(fC1[n-1-i]>=-PointContactR-0.001&&fC1[n-1-i]<=PointContactR+0.001) {
         fV[n-1-i]=V1;
         fIsFixed[n-1-i]=true;
      }
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fV[i]=V0;
      fV[i+n1-1]=V0;
   }
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
   SetupBoundary();
}
