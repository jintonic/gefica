#include "Units.h"
#include "RhoZ.h"
using namespace GeFiCa;

void RhoZ::DoSOR2(int idx)
{
   if (fIsFixed[idx])return; 
   // 2nd-order Successive Over-Relaxation
   double density=fImpurity[idx]*Qe;
   double drm=fdC1m[idx]; // dr_minus
   double drp=fdC1p[idx];
   double dzm=fdC2m[idx];
   double dzp=fdC2p[idx];
   double pzm,pzp,prm,prp; // pzm: potential_z_plus
   if(idx>=fN1)pzm=fV[idx-fN1];
   else pzm=fV[idx+fN1];
   if(idx>=fN-fN1)pzp=fV[idx];
   else pzp=fV[idx+fN1];
   if(idx%fN1==0)prm=fV[idx];
   else prm=fV[idx-1];
   if(idx%fN1==fN1-1)prp=fV[idx];
   else prp=fV[idx+1];
   double tmp=(density/epsilon
         + 1/fC1[idx]*(prp-prm)/(drm+drp) +(prp/drp+prm/drm)*2/(drm+drp)
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
double RhoZ::GetC()
{
   Info("GetC", "Start...");
   CalculatePotential(GeFiCa::kSOR2); // identify undepleted region
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
   CalculatePotential(GeFiCa::kSOR2);
   // set impurity back
   if(fImpurity!=tmpImpurity) delete []fImpurity;
   fImpurity=tmpImpurity;

   // calculate C based on CV^2/2 = epsilon int E^2 dx^3 / 2
   double dV=V0-V1; if(dV<0)dV=-dV;
   double SumofElectricField=0;
   for(int i=0;i<fN;i++) {
      double e1=fE1[i];
      double e2=fE2[i];
      double dr=fdC1p[i];
      double dz=fdC2p[i];
      SumofElectricField+=(e1*e1+e2*e2)*fC1[i]*dr*dz;
   }
   double c=SumofElectricField*2*3.14159*epsilon/dV/dV;
   Info("GetC","%.2f pF",c/pF);
   return c;
}

