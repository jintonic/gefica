#include "Planar1D.h"

using namespace GEFICA;

void Planar1D::SetVoltage(double Voltage)
{
  isbegin[0]=true;
  isbegin[n-1]=true;
  P[n-1]=Voltage;
}
