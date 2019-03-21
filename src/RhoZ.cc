#include "Units.h"
#include "RhoZ.h"
using namespace GeFiCa;

void RhoZ::OverRelaxAt(int idx)
{
   if (fIsFixed[idx])return; 
   // 2nd-order Successive Over-Relaxation
   double density=fImpurity[idx]*Qe;
   double drm=dC1m[idx]; // dr_minus
   double drp=dC1p[idx];
   double dzm=dC2m[idx];
   double dzp=dC2p[idx];
   double pzm,pzp,prm,prp; // pzm: potential_z_plus
   if(idx>=N1)pzm=V[idx-N1];
   else pzm=V[idx+N1];
   if(idx>=fN-N1)pzp=V[idx];
   else pzp=V[idx+N1];
   if(idx%N1==0)prm=V[idx];
   else prm=V[idx-1];
   if(idx%N1==N1-1)prp=V[idx];
   else prp=V[idx+1];
   double tmp=(density/epsilon
         + 1/C1[idx]*(prp-prm)/(drm+drp) +(prp/drp+prm/drm)*2/(drm+drp)
         + (pzp/dzp+pzm/dzm)*2/(dzp+dzm))/
      ((1/drm+1/drp)*2/(drm+drp)+(1/dzp+1/dzm)*2/(dzp+dzm));
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=prm;
   double max=prm;
   if(min>prp)min=prp;
   if (min>pzp)min=pzp;
   if (min>pzm)min=pzm;
   
   //find max
   if(max<prp)max=prp;
   if (max<pzp)max=pzp;
   if (max<pzm)max=pzm;
//if tmp is greater or smaller than max and min, set tmp to it.
   
      //over relax
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
double RhoZ::GetC()
{
   Info("GetC", "Start...");
   SuccessiveOverRelax(); // identify undepleted region
   // set impurity to zero
   double *tmpImpurity=fImpurity;
   for (int i=0;i<fN;i++) {
      if (fImpurity[i]!=0) {
         //impurity not clear,return
         //return -1;
         fImpurity=new double[fN];
         for (int j=0;j<fN;j++) {
            fImpurity[j]=0;
            if (!fIsFixed[j] && !fIsDepleted[j]) fIsFixed[j]=true;
         }
         break;
      }
   }
   // calculate potential without impurity
   SuccessiveOverRelax();
   // set impurity back
   if(fImpurity!=tmpImpurity) delete []fImpurity;
   fImpurity=tmpImpurity;

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   double dV=Bias[0]-Bias[1]; if(dV<0)dV=-dV;
   double SumofElectricField=0;
   for(int i=0;i<fN;i++) {
      double e1=E1[i];
      double e2=E2[i];
      double dr=dC1p[i];
      double dz=dC2p[i];
      SumofElectricField+=(e1*e1+e2*e2)*C1[i]*dr*dz;
   }
   double c=SumofElectricField*2*3.14159*epsilon/dV/dV;
   Info("GetC","%.2f pF",c/pF);
   return c;
}

