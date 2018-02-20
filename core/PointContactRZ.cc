#include "PointContactRZ.h"
#include "iostream"
#include "Units.h"
using namespace GeFiCa;

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
   if (ZLowerBound>=ZUpperBound) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            ZLowerBound, ZUpperBound);
      return;
   }
   double RUpperBound,RLowerBound,PointBegin,PointEnd;
   RUpperBound=Radius;
   RLowerBound=-Radius;
   PointBegin=-PointR;
   PointEnd=PointR;
   double steplength1=(RUpperBound-RLowerBound)/(n1-1);
   double steplength2=(ZUpperBound-ZLowerBound)/(n2-1);
   SetStepLength(steplength1,steplength2);

   // set initial potential values
   for(int i=n;i-->0;) {
      fC1[i]=fC1[i]+RLowerBound;
      fPotential[i]=(V0+V1)/2;
      // set potential for inner electrodes
      if(fC1[i]>PointBegin&&fC1[i]<PointEnd&&fC2[i]<PointDepth) {
         fPotential[i]=V1;
         fIsFixed[i]=true;
      }
   }
   // set potential for outer electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fPotential[i]=V0;
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fPotential[i]=V0;
      fPotential[i+n1-1]=V0;
   }
}
//_____________________________________________________________________________
//
bool PointContactRZ::CalculateField(EMethod method)
{
   if (!fIsLoaded) Initialize();
   return RZ::CalculateField(method);
}
