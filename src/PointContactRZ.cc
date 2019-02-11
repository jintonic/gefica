#include "PointContactRZ.h"
#include "iostream"
#include "Units.h"
using namespace GeFiCa;

void PointContactRZ::BoundaryOnPointcontact()
{
   for(int i=0;i<n;i++)
   {
      //only change fdc1m when it is right close to bound
      if(fC1[i]-PointContactR<fdC1m[i]&&fC1[i]>PointContactR&&fC2[i]<PointContactZ)
      {
         fdC1m[i]=fC1[i]-PointContactR;
      }
      //only change fdc2m when it is right close to bound
      if(fC2[i]-PointContactZ<fdC2m[i]&&fC2[i]>PointContactZ&&fC1[i]<PointContactR)
      {
         fdC2m[i]=fC2[i]-PointContactZ;
      }
   }
}
void PointContactRZ::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (Z0>=Z) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            Z0, Z);
      return;
   }
   double RUpperBound,RLowerBound,PointBegin,PointEnd;
   RUpperBound=Radius;
   PointBegin=0;
   PointEnd=PointContactR;
   double steplength1=(RUpperBound)/(n1-1);
   double steplength2=(Z-Z0)/(n2-1);
   RLowerBound=steplength1/2;//(2.0*Radius/(2*n1-1))/2;
   SetStepLength(steplength1,steplength2);
   for(int i=n;i-->0;) 
   {
      fC1[i]=fC1[i]+RLowerBound;
      fPotential[i]=(V0+V1)/2;
      if(fC1[i]>PointBegin&&fC1[i]<PointEnd&&fC2[i]<PointContactZ)
      {
	    fPotential[i]=V1;
	    fIsFixed[i]=true;
      }
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
      //fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      //fPotential[i]=V0;
      fPotential[i+n1-1]=V0;
   }
   BoundaryOnPointcontact();
}
//_____________________________________________________________________________
//
bool PointContactRZ::CalculatePotential(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return RZ::CalculatePotential(method);
}
