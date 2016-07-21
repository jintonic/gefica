#include "Polar.h"
using namespace GEFICA;

void Polar::CreateGridWithFixedStepLength(double steplength)
{
   XY::CreateGridWithFixedStepLength(steplength);
   for (int i=0;i<n;i++) {
      if(i>n1-1)fC2[i]=fC2[i-n1]*360/n2;
      else fC2[i]=0;
      if(i%n1==0)fC1[i]=0;
      else fC1[i]=fC1[i-1]+steplength;

      fE2[i]=0;
      fDistanceToLeft[i]=360/n2;
      fDistanceToRight[i]=360/n2;
   }
}

double Polar::GetData(double tarx, double tary, int thing)
{
   int idx=FindIdx(tarx,tary,0,n);
   double ab=(tarx-fC1[idx])/fDistanceToNext[idx];
   double aa=1-ab;
   double ba=(tary-fC2[idx])/fDistanceToRight[idx];
   double bb=1-ba;
   double tar0,tar1,tar2,tar3,*tar=NULL;
   switch(thing) {
      case 0:tar= fImpurity;break;
      case 1:tar= fPotential;break;
      case 2:tar= fE1;break;
      case 3:tar= fE2;break;
   }
   tar3=-1;
   tar0=tar[idx];
   if((idx%n1)==n1-1) {
      tar1=tar[idx+1-n1];
      tar3=tar[idx+1];
   } else {
      tar1=tar[idx+1];
   }
   if(idx>n-n1){tar2=0;tar3=0;}
   else {tar2=tar[idx+n1];}
   if (tar3==-1)tar3=tar[idx+n1+1];
   return (tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb;
}
