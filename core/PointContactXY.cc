#include "PointContactXY.h"
using namespace GEFICA;

void PointContactXY::SetVoltage(double dVoltage,double topbegin,double topend)
{
   for(int i=n;i-->n-n1;) {
      fIsFixed[i]=true;
      fPotential[i]=dVoltage;
      if(fC1[n-i]>=topbegin&&fC1[n-i]<=topend) fIsFixed[n-i]=true;
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fPotential[i]=dVoltage;
      fPotential[i+n1-1]=dVoltage;
   }
}
