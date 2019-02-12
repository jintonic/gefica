#include "ReversedCoaxialRZ.h"
#include "iostream"
#include "Units.h"
#include <cmath>
using namespace GeFiCa;
void ReversedCoaxialRZ::Boundary()
{
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
      //right side of hole
      if(fC1[i]-fC2[i]/k1+b1/k1<fdC1m[i]&&fC2[i]>HoleZ)
      {
         fdC1m[i]=fC1[i]-fC2[i]/k1+b1/k1;
      }
      //right of edge
      if(fC1[i]+fC2[i]/k2+b2/k2<fdC1m[i]&&fC2[i]>y4)
      {
         fdC1m[i]=fC1[i]+fC2[i]/k2+b2/k2;
      }
      //left side of hole
      if(-fC1[i]-fC2[i]/k1-b1/k1<fdC1m[i]&&fC2[i]>HoleZ)
      {
         fdC1m[i]=-fC1[i]-fC2[i]/k1-b1/k1;
      }
      //left of edge
      if(-fC1[i]+fC2[i]/k2-b2/k2<fdC1m[i]&&fC2[i]>y4)
      {
         fdC1m[i]=-fC1[i]+fC2[i]/k2-b2/k2;
      }
      //down right side of hole
      if(-fC2[i]+fC1[i]*k1+b1<fdC2m[i]&&fC1[i]>x2)
      {
         fdC2m[i]=-fC1[i]+fC2[i]*k1+b1;
      }
      //down right of edge
      if(-fC2[i]-fC1[i]*k2+b2<fdC2m[i]&&fC1[i]>x2)
      {
         fdC2m[i]=-fC1[i]-fC2[i]*k2+b2;
      }
      //down left side of hole
      if(-fC2[i]-fC1[i]*k1+b1<fdC2m[i]&&fC1[i]>x2)
      {
         fdC2m[i]=-fC1[i]-fC2[i]*k1+b1;
      }
      //down left of edge
      if(-fC2[i]+fC1[i]*k2+b2<fdC2m[i]&&fC1[i]>x2)
      {
         fdC2m[i]=-fC1[i]+fC2[i]*k2+b2;
      }


   }

}
void ReversedCoaxialRZ::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
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
      if(fC1[n-1-i]>=PointContactR-0.001&&fC1[n-1-i]<=-PointContactR+0.001) {
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
     if(((fC2[i]>-k1*(fC1[i])+b1-steplength2/2  && fC2[i]>y2-steplength2/2)||(fC2[i]>-k2*(fC1[i])+b2-steplength2/2))&&fC1[i]<0)
     {
       fIsFixed[i]=true;
       fV[i]=V0;
     }

     if(((fC2[i]>-k1*(fC1[i])+b1+steplength2 && fC2[i]>y2+steplength2)||fC2[i]>-k2*(fC1[i])+b2+steplength2)&&fC1[i]<0)
     {
       fV[i]=0;
     }
     if(((fC2[i]>k1*(fC1[i])+b1-steplength2/2  && fC2[i]>y2-steplength2/2)||(fC2[i]>k2*(fC1[i])+b2-steplength2/2))&&fC1[i]>0)
     {
       fIsFixed[i]=true;
       fV[i]=V0;
     }

     if(((fC2[i]>k1*(fC1[i])+b1+steplength2 && fC2[i]>y2+steplength2)||fC2[i]>k2*(fC1[i])+b2+steplength2)&&fC1[i]>0)
     {
       fV[i]=0;
     }

   }
   //Boundary();
}
//_____________________________________________________________________________
//
bool ReversedCoaxialRZ::CalculatePotential(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return RZ::CalculatePotential(method);
}
