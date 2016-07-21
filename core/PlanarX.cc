#include "PlanarX.h"
using namespace GEFICA;

void PlanarX::SetVoltage(double anode_voltage, double cathode_voltage)
{
   double stepLength=Thickness/n;
   CreateGridWithFixedStepLength(stepLength);
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (cathode_voltage-anode_voltage)/n;
   for (int i=0; i<n; i++) {
      fPotential[i]=anode_voltage+slope*i;
   }
}
