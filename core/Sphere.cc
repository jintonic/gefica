#include "Sphere.h"
using namespace GeFiCa;

void Sphere::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (LowerBound>=UpperBound) {
      Warning("CreateGridWithFixedStepLength",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            LowerBound, UpperBound);
      return;
   }
   double steplength=(UpperBound-LowerBound)/(n1-1);
   SetStepLength(steplength,180/n2,360/(n3));

  for (int i=0;i<n;i=i+n1)
  {
    fPotential[i]=cathode_voltage;
    fPotential[i+n1-1]=anode_voltage;
    fIsFixed[i]=true;
    fIsFixed[i+n1-1]=true;
  }
}
//_____________________________________________________________________________
//
bool Sphere::CalculateField(EMethod method)
{
  if(!fIsLoaded)Initialize();
  return RThetaPhi::CalculateField(method);
}
