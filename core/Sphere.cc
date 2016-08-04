#include "Sphere.h"
using namespace GEFICA;

void Sphere::SetVoltage(double v1,double v2)
{
  for (int i=0;i<n;i=i+n1)
  {
    fPotential[i]=v1;
    fPotential[i+n1-1]=v2;
    fIsFixed[i]=true;
    fIsFixed[i+n1-1]=true;
  }
}
