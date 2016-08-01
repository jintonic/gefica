#include "Polar.h"
using namespace GEFICA;

void Polar::RK2(int idx)
{//need update
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*1.6e12;
   double h2=fDistanceToPrevious[idx];
   double h3=fDistanceToNext[idx];
   double h4=fDistanceToLeft[idx];
   double h1=fDistanceToRight[idx];
   double Pym1,Pyp1,Pxm1,Pxp1;
   if(idx>=n1)Pym1=fPotential[idx-n1];
   else Pym1=fPotential[idx+n-n1];
   if(idx>=n-n1)Pyp1=fPotential[idx-n+n1];
   else Pyp1=fPotential[idx+n1];
   if(idx%n1==0)Pxm1=fPotential[idx];
   else Pxm1=fPotential[idx-1];
   if(idx%n1==n1-1)Pxp1=fPotential[idx];
   else Pxp1=fPotential[idx+1];
   double r=fC1[idx];
   double tmp = (Pxp1/(h3*(h2+h3))+Pxm1/(h2*(h2+h3))
        +Pyp1/r/r/h4/(h1+h4)+Pym1/r/r/h1/(h1+h4)
        -density/epsilon/2+(Pxp1-Pxm1)/2/r/(h2+h3))
       /(1/h3/(h2+h3)+1/h2/(h2+h3)
        +1/r/r/h1/(h1+h4)+1/r/r/h4/(h1+h4));
   fPotential[idx]=Csor*(tmp-fPotential[idx])+fPotential[idx];
   fE1[idx]=(Pxp1-Pxm1)/(h2+h3);
   fE2[idx]=(Pyp1-Pym1)/(h1+h4);
}
void Polar::CreateGridWithFixedStepLength(double steplength)
{
   XY::CreateGridWithFixedStepLength(steplength);
   for (int i=0;i<n;i++) {
      if(i>n1-1)fC2[i]=(double)(int)(i/n1)*2*3.1415926/n2;
      else fC2[i]=0;
      if(i%n1==0)fC1[i]=0;
      else fC1[i]=fC1[i-1]+steplength;

      fE2[i]=0;
      fDistanceToLeft[i]=2*3.1415926/n2;
      fDistanceToRight[i]=2*3.1415926/n2;
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
