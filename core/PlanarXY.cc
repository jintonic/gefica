#include "PlanarXY.h"
using namespace GEFICA;

void PlanarXY::SetVoltage(double dVoltage)
{
   if(isx) {
      for(int i=0;i<x;i++) {
         fIsFixed[i]=true;
         fIsFixed[n-i-1]=true;
         P[i]=dVoltage;
      }
   } else {
      for(int i=0;i<y;i++) {
         fIsFixed[i*x]=true;
         fIsFixed[(i+1)*x-1]=true;
         P[i*x]=dVoltage;
      }
   }
}
