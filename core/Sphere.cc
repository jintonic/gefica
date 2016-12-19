#include "Sphere.h"
using namespace GeFiCa;

void Sphere::Initialize()
{
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
    fPotential[i]=Vpos;
    fPotential[i+n1-1]=Vneg;
    fIsFixed[i]=true;
    fIsFixed[i+n1-1]=true;
  }
}
//_____________________________________________________________________________
//
bool Sphere::CalculateField(EMethod method)
{
  if(!fIsLoaded)Initialize();
  return X::CalculateField(method);
}
