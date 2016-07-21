#include "PlanarX.h"
using namespace GEFICA;

void PlanarX::SetVoltage(double anode_voltage, double cathode_voltage)
{
   double stepLength=Thickness/(n-1);
   CreateGridWithFixedStepLength(stepLength);
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (cathode_voltage-anode_voltage)/n;
   for (int i=0; i<n; i++) {
      fPotential[i]=anode_voltage+slope*i;
   }
}

bool PlanarX::Analyic()
{
   double cnst1=fPotential[0];
   double cnst2=(fPotential[n-1]-fImpurity[n-1]*1.6e-19/2/epsilon*Thickness*Thickness-cnst1)/Thickness;
   for (int i=0; i<n; i++) {
      fPotential[i] = fImpurity[i]*1.6e-19/2/epsilon*fC1[i]*fC1[i]+cnst2*fC1[i]+cnst1;
      fE1[i]=(fPotential[i+1]-fPotential[i-1])/(fDistanceToNext[i]+fDistanceToPrevious[i]);
   }
   return true;
}
