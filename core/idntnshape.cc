#include "idntnshape.h"

using namespace GEFICA;

void idntnshape::SetVoltage(double dVoltage,double r)
{
for(int i=x-1;i<n;i=i+x)
{
  P[i]=dVoltage;
  isbegin[i]=true;
}
for (int i=0;i<n;i=i+x)
{
  isbegin[i]=true;
}
for(int i=0;i<n;i=i+x*y)
{
  int j=i;
  while(C1[j]<r)
  {
    if(j>=x)break;
    isbegin[j]=true;
    j++;
  }
}

}
