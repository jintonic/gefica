#include "Planar2D.h"

using namespace GEFICA;

void Planar2D::SetVoltage(double dVoltage)
{
  if(isx)
  {
    for(int i=0;i<x;i++)
    {
      isbegin[i]=true;
      isbegin[n-i-1]=true;
      P[i]=dVoltage;
    }
  }
  else
  {
    for(int i=0;i<y;i++)
    {
      isbegin[i*x]=true;
      isbegin[(i+1)*x-1]=true;
      P[i*x]=dVoltage;
    }
  }
}
