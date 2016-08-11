#include "PointContactXY.h"
using namespace GEFICA;
//_______________________________
//
//it look like this:
//____________________
//|                  |
//|                  |
//|                  |
//|                  |
//|                  |
//|                  |
//|                  |
//|                  |
//|       ____       |
ClassImp(PointContactXY)
void PointContactXY::initialize()
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
      Warning("CreateGridWithFixedStepLength",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            XLowerBound, XUpperBound);
      return;
   }
  double steplength1=(XUpperBound-XLowerBound)/n1;
  double steplength2=(YUpperBound-YLowerBound)/n2;
  SetStepLength(steplength1,steplength2);
  for(int i=n;i-->n-n1;) {
      fIsFixed[i]=true;
      fPotential[i]=annode_voltage;
      if(fC1[n-i]>=PointBegin&&fC1[n-i]<=PointEnd) 
      {
	fPotential[n-i]=cathode_voltage;
	fIsFixed[n-i]=true;
      }
  }
  for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fPotential[i]=annode_voltage;
      fPotential[i+n1-1]=annode_voltage;
  }
}
bool PointContactXY::CalculateField(EMethod method)
{
  if(!floaded)initialize();
  return X::CalculateField(method);
}
