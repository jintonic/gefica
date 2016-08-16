#include <iostream>
using namespace std;

#include "TrueCoaxial1D.h"
using namespace GEFICA;
//____________________________________________________
//a Planar detector under 1D coordinate system with analyic
ClassImp(TrueCoaxial1D)
void TrueCoaxial1D::initialize()
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
   double steplength=(OuterRadius-InnerRadius)/(n-1);
   SetStepLength(steplength);
   for(int i=n;i-->0;)fC1[i]=fC1[i]+InnerRadius;
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (cathode_voltage-annode_voltage)/(n-1);
   for (int i=0; i<n; i++) fPotential[i]=annode_voltage+slope*i;
}

#include  <cmath>
bool TrueCoaxial1D::Analyic()
{
   double density=fImpurity[1]*1.6e-19;
   double cnst1=(fPotential[n-1]-fPotential[0]-density*(fC1[n-1]*fC1[n-1]-fC1[0]*fC1[0])/epsilon/4)/(log(fC1[n-1]/fC1[0]));
   double cnst2=fPotential[0]-density*fC1[0]*fC1[0]/epsilon/4-cnst1*log(fC1[0]);
   for (int i=0; i<n; i++) {
      fPotential[i] = fImpurity[i]*1.6e-19/4/epsilon*fC1[i]*fC1[i]+cnst1*log(fC1[i])+cnst2;
      fE1[i]=(fPotential[i+1]-fPotential[i-1])/(fDistanceToNext[i]+fDistanceToPrevious[i]);
   }
   return true;
}

bool TrueCoaxial1D::CalculateField(EMethod method)
{
  if(!fIsLoaded)initialize();
  return Rho::CalculateField(method);
}
