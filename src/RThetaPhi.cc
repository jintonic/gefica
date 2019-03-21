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

void RThetaPhi::OverRelaxAt(int idx)
{
   if (fIsFixed[idx])return;

   double density=fImpurity[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double h0=dC3m[idx];
   double h5=dC3p[idx];

   // get potentials of points around point idx
   double pthetam,pthetap,prm,prp,pphip,pphim;
   if (idx<N1*N2) pphim=V[idx+fN-N1*N2];
   else pphim=V[idx-N1*N2];
   if (idx>=fN-N1*N2) pphip=V[idx-(fN-N1*N2)];
   else pphip=V[idx+N1*N2];
   if (idx%(N1*N2)>(N1*N2)-N1-1) {
      if(idx<fN/2) pthetap=V[idx+fN/2];
      else pthetap=V[idx-fN/2];
   } else
      pthetap=V[idx+N1];
   if (idx%(N1*N2)<N1) {
      if(idx<fN/2)pthetam=V[idx+fN/2];
      else pthetam=V[idx-fN/2];
   } else
      pthetam=V[idx-N1];
   if ((idx%(N1*N2))%N1==N1-1) prp=V[idx];
   else prp=V[idx+1];
   if ((idx%(N1*N2))%N1==0) prm=V[idx];
   else prm=V[idx-1];

   double r=C1[idx];
   double O=C2[idx];
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
double RThetaPhi::GetData(double x, double y, double z, double *data)
{
   int idx=FindIdx(x,y,z,0,fN);
   double ab=(x-C1[idx])/dC1p[idx];
   double aa=1-ab;
   double ba=(y-C2[idx])/dC2p[idx];
   double bb=1-ba;
   double ac=(z-C3[idx])/dC3p[idx];
   double ca=1-ac;
   double tar0,tar1,tar2,tar3,tar4,tar5,tar6,tar7;
   if(y==0)return (data[N1*N2-1]+data[N1*N2-1+fN/2])/2;
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=data[idx];
   if(idx>=(fN-N1*N2)) {
      tar4=data[idx-fN+N1*N2];
      tar5=data[idx-fN+N1*N2+1];
      tar6=data[idx-fN+N1*N2+N1];
      tar7=data[idx-fN+N1*N2+N1+1];
   } else
      tar4=data[idx+N1*N2];

   if(idx%(N1*N2)%N1==N1-1) {tar2=0;tar3=0;tar6=0;tar7=0;}
   else {tar2=data[idx+N1];}
   if(idx%(N1*N2)/N1==N2-1) {tar1=0;tar3=0;tar5=0;tar7=0;}
   else {tar1=data[idx+1];}
   if(tar3==-1)tar3=data[idx+N1+1];
   if(tar5==-1)tar5=data[idx+N1*N2+1];
   if(tar6==-1)tar6=data[idx+N1*N2+N1];
   if(tar7==-1)tar7=data[idx+N1*N2+N1+1];

   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac
      +((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
