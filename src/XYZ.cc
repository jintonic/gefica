#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "XYZ.h"
#include "Units.h"
using namespace GeFiCa;

XYZ::XYZ(int nx, int ny, int nz, const char *name, const char *title)
   : XY(nx, ny*nz, name, title)
{ 
   N2=ny; // N2 is set to ny*nz through XY constructor, it is fixed here
   N3=nz;
}
//_____________________________________________________________________________
//
XYZ::~XYZ()
{
   if (E3) delete[] E3;
   if (C3) delete[] C3;
   if (dC3p) delete[] dC3p;
   if (dC3m) delete[] dC3m; 
}
//_____________________________________________________________________________
//
void XYZ::SetStepLength(double steplength1,double steplength2,double steplength3)
{
   XY::SetStepLength(steplength1,steplength2); 
   for (int i=0;i<fN;i++) {
      if(i/(N1*N2)==0) C3[i]=0;
      else C3[i]=C3[i-N1*N2]+steplength3;
      if((i%(N1*N2))/N1!=0)C2[i]=C2[i-N1]+steplength2;
      else C2[i]=0;
      if(i%N1==0)C1[i]=0;
      else C1[i]=C1[i-1]+steplength1;

      E3[i]=0;
      dC3p[i]=steplength3;
      dC3m[i]=steplength3;
   }
}
//_____________________________________________________________________________
//
void XYZ::OverRelaxAt(int idx)
{
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double h0=dC3m[idx];
   double h5=dC3p[idx];
   double pym,pyp,pxm,pxp,pzp,pzm;
   if(idx<N1*N2)pzm=V[idx];
   else pzm=V[idx-N1*N2];
   if(idx>=fN-N1*N2)pzp=V[idx];
   else pzp=V[idx+N1*N2];
   if(idx%(N1*N2)>(N1*N2)-N1-1) pyp=V[idx];
   else pyp=V[idx+N1];
   if(idx%(N1*N2)<N1)pym=V[idx];
   else pym=V[idx-N1];
   if((idx%(N1*N2))%N1==N1-1)pxp=V[idx];
   else pxp=V[idx+1];
   if((idx%(N1*N2))%N1==0)pxm=V[idx];
   else pxm=V[idx-1];

   double tmp= (
         -density/epsilon*h0*h1*h2*h3*h4*h5*(h1+h4)*(h2+h3)*(h0+h5)/2
         +(pxp*h3+pxm*h2)*h0*h1*h4*h5*(h1+h4)*(h0+h5)
         +(pyp*h4+pym*h1)*h0*h2*h3*h5*(h0+h5)*(h2+h3)
         +(pzp*h5+pzm*h0)*h1*h2*h3*h4*(h1+h4)*(h2+h3)	
         )
      /((h0+h5)*(h1+h4)*(h2+h3)*(h0*h1*h4*h5+h0*h2*h3*h5+h1*h2*h3*h4));

   double min=pxm;
   double max=pxm;
   if(min>pxp)min=pxp;
   if (min>pyp)min=pyp;
   if (min>pym)min=pym;
   if (min>pzm)min=pzm;
   if (min>pzm)min=pzm;

   //find max
   if(max<pxp)max=pxp;
   if (max<pyp)min=pyp;
   if (max<pym)max=pym;
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
double XYZ::GetData(double x, double y, double z, double *data)
{
   int idx=FindIdx(x,y,z);
   double ab=(-x+C1[idx])/dC1p[idx];
   double aa=1-ab;
   double ba=(-y+C2[idx])/dC2p[idx];
   double bb=1-ba;
   double ac=(-z+C3[idx])/dC3p[idx];
   double ca=1-ac;
   double tar0,tar1,tar2,tar3,tar4,tar5,tar6,tar7;
   tar3=-1;
   tar5=-1;
   tar6=-1;
   tar7=-1;
   tar0=data[idx];
   if(idx>=(fN-N1*N2)){tar4=0;tar5=0;tar6=0;tar7=0;}
   else{tar4=data[idx-N1*N2];}
   if(idx%(N1*N2)%N1==N1-1){tar2=0;tar3=0;tar6=0;tar7=0;}
   else{tar2=data[idx-N1];}
   if(idx%(N1*N2)/N1==N2-1){tar1=0;tar3=0;tar5=0;tar7=0;}
   else{tar1=data[idx+1];}
   if(tar3==-1)tar3=data[idx-N1-1];
   if(tar5==-1)tar5=data[idx-N1*N2-1];
   if(tar6==-1)tar6=data[idx-N1*N2-N1];
   if(tar7==-1)tar7=data[idx-N1*N2-N1-1];
   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac
      +((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
//_____________________________________________________________________________
//
bool XYZ::CalculateField(int idx)
{
   if (!XY::CalculateField(idx)) return false;
   if (dC3p[idx]==0 || dC3m[idx]==0) return false;

   if (idx<N1*N2) // C3 lower border
      E3[idx]=(V[idx]-V[idx+N1])/dC3p[idx];
   else if (idx>=fN-N1*N2) // C3 upper border
      E3[idx]=(V[idx-N1]-V[idx])/dC3m[idx];
   else { // bulk
      E3[idx]=(V[idx-N1]-V[idx+N1])/(dC3m[idx]+dC3p[idx]);
   }
   return true;
}
