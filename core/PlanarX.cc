#include "PlanarX.h"
using namespace GEFICA;

void PlanarX::SetVoltage(double Voltage)
{
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   P[n-1]=Voltage;
}
