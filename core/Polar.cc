#include "Polar.h"
using namespace GEFICA;
void Polar::Create(double steplength)
{
  Field2D::Create(steplength);
  for (int i=0;i<n;i++)
  {
    if(i>x-1)C2[i]=C2[i-x]*360/y;
    else C2[i]=0;
    if(i%x==0)C1[i]=0;
    else C1[i]=C1[i-1]+steplength;

    E2[i]=0;
    StepLeft[i]=360/y;
    StepRight[i]=360/y;
  }
}
double Polar::GetData(double tarx, double tary, int thing)
{
  int idx=FindIdx(tarx,tary,0,n);
  double ab=(tarx-C1[idx])/StepNext[idx];
  double aa=1-ab;
  double ba=(tary-C2[idx])/StepRight[idx];
  double bb=1-ba;
  double tar0,tar1,tar2,tar3,*tar=NULL;
  switch(thing)
  {
    case 0:tar= Impurity;break;
    case 1:tar= P;break;
    case 2:tar= E1;break;
    case 3:tar= E2;break;
  }
  tar3=-1;
  tar0=tar[idx];
  if((idx%x)==x-1)
  {
    tar1=tar[idx+1-x];
    tar3=tar[idx+1];
  }
  else {tar1=tar[idx+1];}
  if(idx>n-x){tar2=0;tar3=0;}
  else {tar2=tar[idx+x];}
  if (tar3==-1)tar3=tar[idx+x+1];
  return (tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb;
}
