#include <TF3.h>
#include <TFile.h>
#include <TChain.h>
#include <TVectorD.h>

#include "XYZ.h"
#include "Units.h"
using namespace GeFiCa;

XYZ::XYZ(int nx, int ny, int nz, const char *name, const char *title)
   : XY(nx, ny*nz, name, title), n3(nz)
{ 
   n2=ny; // n2 is set to ny*nz through XY constructor, it is fixed here
}
//_____________________________________________________________________________
//
XYZ::~XYZ()
{
   if (fE3) delete[] fE3;
   if (fC3) delete[] fC3;
   if (fdC3p) delete[] fdC3p;
   if (fdC3m) delete[] fdC3m; 
}
//_____________________________________________________________________________
//
void XYZ::SetStepLength(double steplength1,double steplength2,double steplength3)
{
   XY::SetStepLength(steplength1,steplength2); 
   for (int i=0;i<n;i++) {
      if(i/(n1*n2)==0) fC3[i]=0;
      else fC3[i]=fC3[i-n1*n2]+steplength3;
      if((i%(n1*n2))/n1!=0)fC2[i]=fC2[i-n1]+steplength2;
      else fC2[i]=0;
      if(i%n1==0)fC1[i]=0;
      else fC1[i]=fC1[i-1]+steplength1;

      fE3[i]=0;
      fdC3p[i]=steplength3;
      fdC3m[i]=steplength3;
   }
}
//_____________________________________________________________________________
//
void XYZ::DoSOR2(int idx)
{
   if (fIsFixed[idx])return;
   double density=-fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h4=fdC2m[idx];
   double h1=fdC2p[idx];
   double h0=fdC3m[idx];
   double h5=fdC3p[idx];
   double pym,pyp,pxm,pxp,pzp,pzm;
   if(idx<n1*n2)pzm=fV[idx];
   else pzm=fV[idx-n1*n2];
   if(idx>=n-n1*n2)pzp=fV[idx];
   else pzp=fV[idx+n1*n2];
   if(idx%(n1*n2)>(n1*n2)-n1-1) pyp=fV[idx];
   else pyp=fV[idx+n1];
   if(idx%(n1*n2)<n1)pym=fV[idx];
   else pym=fV[idx-n1];
   if((idx%(n1*n2))%n1==n1-1)pxp=fV[idx];
   else pxp=fV[idx+1];
   if((idx%(n1*n2))%n1==0)pxm=fV[idx];
   else pxm=fV[idx-1];

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
int XYZ::FindIdx(double tarx, double tary ,double tarz,int begin,int end)
{
   //search using binary search
   if(begin>=end)return XY::FindIdx(tarx,tary,begin,begin+n1*n2-1);
   int mid=((begin/(n1*n2)+end/(n1*n2))/2)*n1*n2;
   if(fC3[mid]>=tarz)return FindIdx(tarx,tary,tarz,begin,mid);
   else return FindIdx(tarx,tary,tarz,mid+1,end);
}
//_____________________________________________________________________________
//
double XYZ::GetData(double tarx, double tary, double tarz, EOutput output )
{
   //get item with number: 0:Impurity 1:Potential 2: Ex 3:Ey 4:Ez

   int idx=FindIdx(tarx,tary,tarz,0,n);
   double ab=(-tarx+fC1[idx])/fdC1p[idx];
   double aa=1-ab;
   double ba=(-tary+fC2[idx])/fdC2p[idx];
   double bb=1-ba;
   double ac=(-tarz+fC3[idx])/fdC3p[idx];
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
   if(idx>=(n-n1*n2)){tar4=0;tar5=0;tar6=0;tar7=0;}
   else{tar4=tar[idx-n1*n2];}
   if(idx%(n1*n2)%n1==n1-1){tar2=0;tar3=0;tar6=0;tar7=0;}
   else{tar2=tar[idx-n1];}
   if(idx%(n1*n2)/n1==n2-1){tar1=0;tar3=0;tar5=0;tar7=0;}
   else{tar1=tar[idx+1];}
   if(tar3==-1)tar3=tar[idx-n1-1];
   if(tar5==-1)tar5=tar[idx-n1*n2-1];
   if(tar6==-1)tar6=tar[idx-n1*n2-n1];
   if(tar7==-1)tar7=tar[idx-n1*n2-n1-1];
   return ((tar0*aa+tar1*ab)*ba+(tar2*aa+tar3*ab)*bb)*ac
      +((tar4*aa+tar5*ab)*ba+(tar6*aa+tar7*ab)*bb)*ca;
}
//_____________________________________________________________________________
//
bool XYZ::CalculateField(int idx)
{
   if (!XY::CalculateField(idx)) return false;
   if (fdC3p[idx]==0 || fdC3m[idx]==0) return false;

   if (idx<n1*n2) // C3 lower border
      fE3[idx]=(fV[idx]-fV[idx+n1])/fdC3p[idx];
   else if (idx>=n-n1*n2) // C3 upper border
      fE3[idx]=(fV[idx]-fV[idx-n1])/fdC3m[idx];
   else { // bulk
      fE3[idx]=(fV[idx-n1]-fV[idx+n1])/(fdC3m[idx]+fdC3p[idx]);
   }
   return true;
}
