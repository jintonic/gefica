#include "Units.h"
#include "RhoPhiZ.h"
using namespace GeFiCa;

void RhoPhiZ::OverRelaxAt(int idx)
{//need update
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double h0=dC3m[idx];
   double h5=dC3p[idx];
   double pphim,pphip,prhom,prhop,pzp,pzm;
   if(idx<N1*N2)pzm=V[idx];
   else pzm=V[idx-N1*N2];
   if(idx>=fN-N1*N2)pzp=V[idx];
   else pzp=V[idx+N1*N2];
   if(idx%(N1*N2)>(N1*N2)-N1-1) pphip=V[idx-N1*N2+N1];
   else pphip=V[idx+N1];
   if(idx%(N1*N2)<N1)pphim=V[idx+N1*N2-N1];
   else pphim=V[idx-N1];
   if((idx%(N1*N2))%N1==N1-1)prhop=V[idx];
   else prhop=V[idx+1];
   if((idx%(N1*N2))%N1==0)prhom=V[idx];
   else prhom=V[idx-1];
   double r=C1[idx];
   double tmp= (-density/2/epsilon+(prhop-prhom)/(2*r*(h2+h3))+prhop/h3/(h2+h3)
         +prhom/h2/(h2+h3)+pphip/h4/(h1+h4)/r/r
         +pphim/h1/(h1+h4)/r/r+pzp/h5/(h0+h5)+pzm/h0/(h0+h5))
      /(1/h2/(h2+h3)+1/h3/(h2+h3)+1/h4/(h1+h4)/r/r
            +1/h1/(h1+h4)/r/r+1/h0/(h0+h5)+1/h5/(h0+h5));
   V[idx]=RelaxationFactor*(tmp-V[idx])+V[idx];
   double min=prhom;
   double max=prhom;
   if(min>prhop)min=prhop;
   if (min>pphip)min=pphip;
   if (min>pphim)min=pphim;
   if (min>pzm)min=pzm;
   if (min>pzm)min=pzm;

   //find max
   if(max<prhop)max=prhop;
   if (max<pphip)min=pphip;
   if (max<pphim)max=pphim;
   if (max<pzm)max=pzm;
   if (max<pzm)max=pzm;
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
double RhoPhiZ::GetData(double x, double y, double z, double *data)
{
   //0:Impurity 1:Potential 2:E1 3:E2 3:E3
   int idx=FindIdx(x,y,z,0,fN);
   double ab=(x-C1[idx])/dC1p[idx];
   double aa=1-ab;
   double ba=(y-C2[idx])/dC2p[idx];
   double bb=1-ba;
   double ac=(z-C3[idx])/dC3p[idx];
   double ca=1-ac;
   double tar0,tar1,tar2,tar3,tar4,tar5,tar6,tar7;
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=data[idx];
   if(idx>=(fN-N1*N2)){tar4=0;tar5=0;tar6=0;tar7=0;}
   else{tar4=data[idx+N1*N2];}
   if(idx%(N1*N2)%N1==N1-1){tar2=0;tar3=0;tar6=0;tar7=0;}
   else{tar2=data[idx+N1];}
   if(idx%(N1*N2)/N1==N2-1){tar1=0;tar3=0;tar5=0;tar7=0;}
   else{tar1=data[idx+1];}
   if(tar3==-1)tar3=data[idx+N1+1];
   if(tar5==-1)tar5=data[idx+N1*N2+1];
   if(tar6==-1)tar6=data[idx+N1*N2+N1];
   if(tar7==-1)tar7=data[idx+N1*N2+N1+1];
   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac
      +((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
