#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "RhoPhiZ.h"
#include "Units.h"
using namespace GeFiCa;

void RhoPhiZ::DoSOR2(int idx)
{//need update
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h4=fdC2m[idx];
   double h1=fdC2p[idx];
   double h0=fdC3m[idx];
   double h5=fdC3p[idx];
   double pphim,pphip,prhom,prhop,pzp,pzm;
   if(idx<fN1*fN2)pzm=fV[idx];
   else pzm=fV[idx-fN1*fN2];
   if(idx>=fN-fN1*fN2)pzp=fV[idx];
   else pzp=fV[idx+fN1*fN2];
   if(idx%(fN1*fN2)>(fN1*fN2)-fN1-1) pphip=fV[idx-fN1*fN2+fN1];
   else pphip=fV[idx+fN1];
   if(idx%(fN1*fN2)<fN1)pphim=fV[idx+fN1*fN2-fN1];
   else pphim=fV[idx-fN1];
   if((idx%(fN1*fN2))%fN1==fN1-1)prhop=fV[idx];
   else prhop=fV[idx+1];
   if((idx%(fN1*fN2))%fN1==0)prhom=fV[idx];
   else prhom=fV[idx-1];
   double r=fC1[idx];
   double tmp= (-density/2/epsilon+(prhop-prhom)/(2*r*(h2+h3))+prhop/h3/(h2+h3)
         +prhom/h2/(h2+h3)+pphip/h4/(h1+h4)/r/r
         +pphim/h1/(h1+h4)/r/r+pzp/h5/(h0+h5)+pzm/h0/(h0+h5))
      /(1/h2/(h2+h3)+1/h3/(h2+h3)+1/h4/(h1+h4)/r/r
            +1/h1/(h1+h4)/r/r+1/h0/(h0+h5)+1/h5/(h0+h5));
   fV[idx]=Csor*(tmp-fV[idx])+fV[idx];
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
double RhoPhiZ::GetData(double tarx, double tary, double tarz, EOutput output)
{
   //0:Impurity 1:Potential 2:E1 3:E2 3:E3
   int idx=FindIdx(tarx,tary,tarz,0,fN);
   double ab=(tarx-fC1[idx])/fdC1p[idx];
   double aa=1-ab;
   double ba=(tary-fC2[idx])/fdC2p[idx];
   double bb=1-ba;
   double ac=(tarz-fC3[idx])/fdC3p[idx];
   double ca=1-ac;
   double tar0,tar1,tar2,tar3,tar4,tar5,tar6,tar7,*tar=NULL;
   switch(output)
   {
      case 0:tar= fImpurity;break;
      case 1:tar= fV;break;
      case 2:tar= fE1;break;
      case 3:tar= fE2;break;
      case 4:tar= fE3;break;
   }
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=tar[idx];
   if(idx>=(fN-fN1*fN2)){tar4=0;tar5=0;tar6=0;tar7=0;}
   else{tar4=tar[idx+fN1*fN2];}
   if(idx%(fN1*fN2)%fN1==fN1-1){tar2=0;tar3=0;tar6=0;tar7=0;}
   else{tar2=tar[idx+fN1];}
   if(idx%(fN1*fN2)/fN1==fN2-1){tar1=0;tar3=0;tar5=0;tar7=0;}
   else{tar1=tar[idx+1];}
   if(tar3==-1)tar3=tar[idx+fN1+1];
   if(tar5==-1)tar5=tar[idx+fN1*fN2+1];
   if(tar6==-1)tar6=tar[idx+fN1*fN2+fN1];
   if(tar7==-1)tar7=tar[idx+fN1*fN2+fN1+1];
   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac+((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
