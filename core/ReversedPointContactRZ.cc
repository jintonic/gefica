#include "ReversedPointContactRZ.h"
#include "iostream"
#include "Units.h"
#include <cmath>
using namespace GeFiCa;

void ReversedPointContactRZ::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (RLowerBound>=RUpperBound||ZLowerBound>=ZUpperBound) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            RLowerBound, RUpperBound);
      return;
   }
   double steplength1=(RUpperBound-RLowerBound)/(n1-1);
   double steplength2=(ZUpperBound-ZLowerBound)/(n2-1);
   std::cout<<steplength1<<std::endl; 
   SetStepLength(steplength1,steplength2);
   double x1=OutterRadiusHole,
	  y1=ZUpperBound,
	  x2=InnerRadiusHole,
	  y2=ZUpperBound-DHole,
	  x3=RUpperBound-removedConnorradius,
	  y3=ZUpperBound,
	  x4=RUpperBound,
	  y4=ZUpperBound-removedConnorheight;
   double k1=(y1-y2)/(x1-x2);
   double b1=y1-k1*x1;
   double k2=(y3-y4)/(x3-x4);
   double b2=(y3-k2*x3);

   for(int i=n;i-->0;) 
   {
      fC1[i]=fC1[i]+RLowerBound;
      fPotential[i]=(V0+V1)/2;
   }
   // set potential for electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fPotential[i]=V0;
      if(fC1[n-1-i]>=PointBegin-0.001&&fC1[n-1-i]<=PointEnd+0.001) {
         fPotential[n-1-i]=V1;
         fIsFixed[n-1-i]=true;
      }
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fPotential[i]=V0;
      fPotential[i+n1-1]=V0;
   }
   for (int i=0;i<n;i++)
   {
     if(((fC2[i]>-k1*(fC1[i])+b1-steplength2/2  && fC2[i]>y2-steplength2/2)||(fC2[i]>-k2*(fC1[i])+b2-steplength2/2))&&fC1[i]<0)
     {
       fIsFixed[i]=true;
       fPotential[i]=V0;
     }

     if(((fC2[i]>-k1*(fC1[i])+b1+steplength2 && fC2[i]>y2+steplength2)||fC2[i]>-k2*(fC1[i])+b2+steplength2)&&fC1[i]<0)
     {
       fPotential[i]=0;
     }
     if(((fC2[i]>k1*(fC1[i])+b1-steplength2/2  && fC2[i]>y2-steplength2/2)||(fC2[i]>k2*(fC1[i])+b2-steplength2/2))&&fC1[i]>0)
     {
       fIsFixed[i]=true;
       fPotential[i]=V0;
     }

     if(((fC2[i]>k1*(fC1[i])+b1+steplength2 && fC2[i]>y2+steplength2)||fC2[i]>k2*(fC1[i])+b2+steplength2)&&fC1[i]>0)
     {
       fPotential[i]=0;
     }

   }
}
//_____________________________________________________________________________
//
bool ReversedPointContactRZ::CalculateField(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return RZ::CalculateField(method);
}
