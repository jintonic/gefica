#include "RhoPhi.h"
#include "Units.h"
using namespace GeFiCa;
#include <iostream>
using namespace std;
void RhoPhi::OverRelaxAt(int idx)
{
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double r=C1[idx];
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double pphim,pphip,prhom,prhop;
   if(idx>=N1)pphim=V[idx-N1];
   else pphim=V[idx+fN-N1];
   if(idx>=fN-N1)pphip=V[idx-fN+N1];
   else pphip=V[idx+N1];
   if(idx%N1==0)prhom=V[idx];
   else prhom=V[idx-1];
   if(idx%N1==N1-1)prhop=V[idx];
   else prhop=V[idx+1];
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
   //V[idx]=RelaxationFactor*(tmp-V[idx])+V[idx];
   //if need calculate depleted voltage
   double oldP=V[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;
   if(tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if(tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||Bias[0]==Bias[1]) V[idx]=tmp;
}
//_____________________________________________________________________________
//
double RhoPhi::GetData(double x, double y, double z, double *data)
{
   //0:Impurity 1:Potential 2:E1 3:E2
   int idx=FindIdx(x,y,0,fN);
   double ab=(x-C1[idx])/dC1p[idx];
   double aa=1-ab;
   double ba=(y-C2[idx])/dC2p[idx];
   double bb=1-ba;
   double tar0,tar1,tar2,tar3;
   tar3=-1;
   tar0=data[idx];
   if((idx%N1)==N1-1) {
      tar1=data[idx+1-N1];
      tar3=data[idx+1];
   } else {
      tar1=data[idx+1];
   }
   if(idx>fN-N1){tar2=0;tar3=0;}
   else {tar2=data[idx+N1];}
   if (tar3==-1)tar3=data[idx+N1+1];
   return (tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb;
}
