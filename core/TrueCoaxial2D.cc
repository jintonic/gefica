#include "TrueCoaxial2D.h"
using namespace GEFICA;
//____________________________________________________
//a True Coaxial detector under 2D coordinate system
ClassImp(TrueCoaxial2D)
void TrueCoaxial2D::initialize()
{
  // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(OuterRadius-InnerRadius)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (InnerRadius>=OuterRadius) {
      Warning("CreateGridWithFixedStepLength",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            InnerRadius, OuterRadius);
      return;
   }
   double steplength=(OuterRadius-InnerRadius)/(n1-1);
   SetStepLength(steplength,3.14159265/n2);

  for (int i=0;i<n;i=i+n1)
  {
    fPotential[i]=cathode_voltage;
    fPotential[i+n1-1]=anode_voltage;
    fIsFixed[i]=true;
    fIsFixed[i+n1-1]=true;
  }
}
bool TrueCoaxial2D::CalculateField(EMethod method)
{
  if(!fIsLoaded)initialize();
  return RhoPhi::CalculateField(method);
}
