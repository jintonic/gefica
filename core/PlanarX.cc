#include "PlanarX.h"
using namespace GEFICA;

void PlanarX::SetVoltage(double anode_voltage, double cathode_voltage)
{
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   fPotential[n-1]=(anode_voltage-cathode_voltage);
}
