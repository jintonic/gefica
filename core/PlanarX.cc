#include "PlanarX.h"
using namespace GEFICA;

void PlanarX::SetVoltage(double Voltage)
{
   isbegin[0]=true;
   isbegin[n-1]=true;
   P[n-1]=Voltage;
}