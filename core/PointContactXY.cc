#include "PointContactXY.h"
using namespace GeFiCa;

void PointContactXY::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (XLowerBound>=XUpperBound||YLowerBound>=YUpperBound) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            XLowerBound, XUpperBound);
      return;
   }
   double steplength1=(XUpperBound-XLowerBound)/n1;
   double steplength2=(YUpperBound-YLowerBound)/n2;
   SetStepLength(steplength1,steplength2);
   // set potential for electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fPotential[i]=annode_voltage;
      if(fC1[n-1-i]>=PointBegin&&fC1[n-1-i]<=PointEnd) {
         fPotential[n-1-i]=cathode_voltage;
         fIsFixed[n-1-i]=true;
      }
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fPotential[i]=annode_voltage;
      fPotential[i+n1-1]=annode_voltage;
   }
}
//_____________________________________________________________________________
//
bool PointContactXY::CalculateField(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return XY::CalculateField(method);
}
