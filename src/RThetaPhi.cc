#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "RThetaPhi.h"
#include "Units.h"
using namespace GeFiCa;

#include <cmath>
#include <iostream>
using namespace std;

void RThetaPhi::DoSOR2(int idx)
{
   if (fIsFixed[idx])return;

   double density=fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h4=fdC2m[idx];
   double h1=fdC2p[idx];
   double h0=fdC3m[idx];
   double h5=fdC3p[idx];

   // get potentials of points around point idx
   double pthetam,pthetap,prm,prp,pphip,pphim;
   if (idx<fN1*fN2) pphim=fV[idx+fN-fN1*fN2];
   else pphim=fV[idx-fN1*fN2];
   if (idx>=fN-fN1*fN2) pphip=fV[idx-(fN-fN1*fN2)];
   else pphip=fV[idx+fN1*fN2];
   if (idx%(fN1*fN2)>(fN1*fN2)-fN1-1) {
      if(idx<fN/2) pthetap=fV[idx+fN/2];
      else pthetap=fV[idx-fN/2];
   } else
      pthetap=fV[idx+fN1];
   if (idx%(fN1*fN2)<fN1) {
      if(idx<fN/2)pthetam=fV[idx+fN/2];
      else pthetam=fV[idx-fN/2];
   } else
      pthetam=fV[idx-fN1];
   if ((idx%(fN1*fN2))%fN1==fN1-1) prp=fV[idx];
   else prp=fV[idx+1];
   if ((idx%(fN1*fN2))%fN1==0) prm=fV[idx];
   else prm=fV[idx-1];

   double r=fC1[idx];
   double O=fC2[idx];
   double tmp = (density/epsilon/2
         +(prp-prm)/r/(h2+h3)
         +(pthetap-pthetam)/r/r/(h1+h4)/sin(O)*cos(O)/2
         +(prp/h3+prm/h2)/(h3+h2)
         +(pthetap/h1+pthetam/h4)/r/r/(h4+h1)
         +(pphip/h5+pphim/h0)/r/r/sin(O)/sin(O)/(h0+h5))
      /(     1/h2/h3
            +1/r/r/h1/h4
            +1/r/r/sin(O)/sin(O)/h5/h0);
   double min=prm;
   double max=prm;
   if(min>prp)min=prp;
   if (min>pphip)min=pphip;
   if (min>pphim)min=pphim;
   if (min>pthetam)min=pthetam;
   if (min>pthetam)min=pthetam;

   //find max
   if(max<prp)max=prp;
   if (max<pphip)max=pphip;
   if (max<pphim)max=pphim;
   if (max<pthetam)max=pthetam;
   if (max<pthetam)max=pthetam;
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
double RThetaPhi::GetData(double tarx, double tary, double tarz, EOutput output)
{
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
   if(tary==0)return (tar[fN1*fN2-1]+tar[fN1*fN2-1+fN/2])/2;
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=tar[idx];
   if(idx>=(fN-fN1*fN2)) {
      tar4=tar[idx-fN+fN1*fN2];
      tar5=tar[idx-fN+fN1*fN2+1];
      tar6=tar[idx-fN+fN1*fN2+fN1];
      tar7=tar[idx-fN+fN1*fN2+fN1+1];
   } else
      tar4=tar[idx+fN1*fN2];

   if(idx%(fN1*fN2)%fN1==fN1-1) {tar2=0;tar3=0;tar6=0;tar7=0;}
   else {tar2=tar[idx+fN1];}
   if(idx%(fN1*fN2)/fN1==fN2-1) {tar1=0;tar3=0;tar5=0;tar7=0;}
   else {tar1=tar[idx+1];}
   if(tar3==-1)tar3=tar[idx+fN1+1];
   if(tar5==-1)tar5=tar[idx+fN1*fN2+1];
   if(tar6==-1)tar6=tar[idx+fN1*fN2+fN1];
   if(tar7==-1)tar7=tar[idx+fN1*fN2+fN1+1];

   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac
      +((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
