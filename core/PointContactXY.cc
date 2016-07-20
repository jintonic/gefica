#include "PointContactXY.h"
using namespace GEFICA;

void PointContactXY::SetVoltage(double dVoltage,double topbegin,double topend)
{
   for(int i=n;i-->n-x;) {
      fIsFixed[i]=true;
      P[i]=dVoltage;
      if(C1[n-i]>=topbegin&&C1[n-i]<=topend) fIsFixed[n-i]=true;
   }
   for(int i=0;i<n-x;i=i+x) {
      fIsFixed[i]=true;
      fIsFixed[i+x-1]=true;
      P[i]=dVoltage;
      P[i+x-1]=dVoltage;
   }
}
