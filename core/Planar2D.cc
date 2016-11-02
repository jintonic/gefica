#include "Planar2D.h"
using namespace GeFiCa;

void Planar2D::Initialize()
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
  XY::SetStepLength(steplength1,steplength2);
   for(int i=0;i<n2;i++) {
      fIsFixed[i*n1]=true;
      fIsFixed[(i+1)*n1-1]=true;
      fPotential[i*n1]=cathode_voltage;
      fPotential[(i+1)*n1-1]=annode_voltage;
   }
}
//_____________________________________________________________________________
//
bool Planar2D::CalculateField(EMethod method)
{
  if(!fIsLoaded)Initialize();
  return XY::CalculateField(method);
}
