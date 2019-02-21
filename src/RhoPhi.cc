#include "RhoPhi.h"
#include "Units.h"
using namespace GeFiCa;
#include <iostream>
using namespace std;
void RhoPhi::DoSOR2(int idx)
{
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double r=fC1[idx];
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h4=fdC2m[idx];
   double h1=fdC2p[idx];
   double pphim,pphip,prhom,prhop;
   if(idx>=n1)pphim=fV[idx-n1];
   else pphim=fV[idx+n-n1];
   if(idx>=n-n1)pphip=fV[idx-n+n1];
   else pphip=fV[idx+n1];
   if(idx%n1==0)prhom=fV[idx];
   else prhom=fV[idx-1];
   if(idx%n1==n1-1)prhop=fV[idx];
   else prhop=fV[idx+1];
   double tmp = (prhop/(h3*(h2+h3))+prhom/(h2*(h2+h3))
         +pphip/r/r/h4/(h1+h4)+pphim/r/r/h1/(h1+h4)
         -density/epsilon/2+(prhop-prhom)/2/r/(h2+h3))
      /(1/h3/(h2+h3)+1/h2/(h2+h3)
            +1/r/r/h1/(h1+h4)+1/r/r/h4/(h1+h4));
   double min=prhom;
   double max=prhom;
   if(min>prhop)min=prhop;
   if (min>pphip)min=pphip;
   if (min>pphim)min=pphim;

   //find max
   if(max<prhop)max=prhop;
   if (max<pphip)max=pphip;
   if (max<pphim)max=pphim;
   //if tmp is greater or smaller than max and min, set tmp to it.
   //fV[idx]=Csor*(tmp-fV[idx])+fV[idx];
   //if need calculate depleted voltage
   double oldP=fV[idx];
   tmp=Csor*(tmp-oldP)+oldP;
   if(tmp<min) {
      fV[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      fV[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||V0==V1) fV[idx]=tmp;
}
//_____________________________________________________________________________
//
double RhoPhi::GetData(double tarx, double tary, EOutput output)
{
   //0:Impurity 1:Potential 2:E1 3:E2
   int idx=FindIdx(tarx,tary,0,n);
   double ab=(tarx-fC1[idx])/fdC1p[idx];
   double aa=1-ab;
   double ba=(tary-fC2[idx])/fdC2p[idx];
   double bb=1-ba;
   double tar0,tar1,tar2,tar3,*tar=NULL;
   switch(output) {
      case 0:tar= fImpurity;break;
      case 1:tar= fV;break;
      case 2:tar= fE1;break;
      case 3:tar= fE2;break;
      default:break;
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
