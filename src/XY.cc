#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "XY.h"
#include "Units.h"
using namespace GeFiCa;

XY::XY(int nx, int ny, const char *name, const char *title)
   : X(nx*ny, name, title)
{
   fN1=nx; // fN1 is set to nx*ny through X constructor, it is fixed here
   fN2=ny;
   fEgraphs = new TMultiGraph("gEs","electric field lines");
}
//_____________________________________________________________________________
//
XY::~XY()
{
   if (fE2) delete[] fE2;
   if (fC2) delete[] fC2;
   if (fdC2m) delete[] fdC2m;
   if (fdC2p) delete[] fdC2p;
   if (fEgraphs) delete fEgraphs;
}
//_____________________________________________________________________________
//
void XY::SetStepLength(double steplength1,double steplength2)
{
   //set field step length
   X::SetStepLength(steplength1);
   for (int i=0;i<fN;i++) {
      if(i>fN1-1)fC2[i]=fC2[i-fN1]+steplength2;
      else fC2[i]=0;
      if(i%fN1==0)fC1[i]=0;
      else fC1[i]=fC1[i-1]+steplength1;

      fE2[i]=0;
      fdC2m[i]=steplength2;
      fdC2p[i]=steplength2;
   }
}
//_____________________________________________________________________________
//
#include <iostream>
using namespace std;
void XY::DoSOR2(int idx)
{
   // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*Qe;
   double h2=fdC1m[idx];
   double h3=fdC1p[idx];
   double h4=fdC2m[idx];
   double h1=fdC2p[idx];
   double pym,pyp,pxm,pxp;
   if(idx>=fN1)pym=fV[idx-fN1];
   else pym=fV[idx];
   if(idx>=fN-fN1)pyp=fV[idx];
   else pyp=fV[idx+fN1];
   if(idx%fN1==0)pxm=fV[idx];
   else pxm=fV[idx-1];
   if(idx%fN1==fN1-1)pxp=fV[idx];
   else pxp=fV[idx+1];
   double tmp=(density/epsilon+(pxp/h2+pxm/h3)*2/(h2+h3)+(pyp/h1+pym/h4)*2/(h1+h4))/
      ((1/h2+1/h3)*2/(h2+h3)+(1/h1+1/h4)*2/(h1+h4));
   //find minmium and maxnium of all five grid, the new one should not go overthem.
   //find min
   double min=pxm;
   double max=pxm;
   if(min>pxp)min=pxp;
   if (min>pyp)min=pyp;
   if (min>pym)min=pym;
   
   //find max
   if(max<pxp)max=pxp;
   if (max<pyp)max=pyp;
   if (max<pym)max=pym;
   //if tmp is greater or smaller than max and min, set tmp to it.
   
   //fV[idx]=RelaxationFactor*(tmp-fV[idx])+fV[idx];
   //if need calculate depleted voltage
   double oldP=fV[idx];
   tmp=RelaxationFactor*(tmp-oldP)+oldP;
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
double XY::GetData(double x, double y, double z, double *data)
{
   // +-aa--+--ab--+(fC1[idx], fC2[idx])
   // |     ^      |
   // |  dyp|      ba
   // |     |(x,y) |
   // +<----+----->+
   // | dxm | dxp  |
   // |     |      bb
   // |     v dym  |
   // +-----+------+
   int idx=FindIdx(x,y);

   bool tr=true; // value of top right grid point 
   bool tl=false; // value of top left grid point
   bool br=false; // value of bottom right grid point
   bool bl=false; // value of bottom left grid point
   if (idx%fN1==0) {  } // left boundary
   else tl=true;
   if (idx<fN1) {  } // bottom boundary
   else br=true;
   if (tlv&&b) bl=true; // neither left nor bottom boundary

   double topmiddlevalue,bottommdlevalue,bottommiddiledistencetoabove,topmiddledistencetobelow;
   if(tr&&tl )
   {
      double trv,tlv,aa,ab;
      double boundaryoint=fdC1m[idx]<fdC1p[idx-1] ? fC1[idx]-fdC1m[idx] : fC1[idx-1]+fdC1p[idx-1];
      if(x<boundaryoint)
      {
         tlv=fV[idx-1];
         if(fIsFixed[idx-1])
         {
            trv=fV[idx-1];
         }
         else//right fix or none fix
         {
            trv=fV[idx];
         }
         aa=x-fC1[idx-1];
         ab=boundaryoint-x;
      }
      else
      {
         trv=fV[idx];
         if(fIsFixed[idx])
         {
            tlv=fV[idx];
         }
         else
         {
            tlv=fV[idx-1];
         }
      }

   }
   else topmiddlevalue=fV[idx];
   if(br&&bl)
   {
   }
   else 
        
   double ab=(-x+fC1[idx])/fdC1m[idx];
   double aa=1-ab;
   double ba=(-y+fC2[idx])/fdC2m[idx];
   double bb=1-ba;

   if (gDebug>0) {
      Info("GetData","fC1[%d]=%f, ",idx,fC1[idx]);
   }
   return (tlv*ab+trv*aa)*bb+(blv*ab+brv*aa)*ba;
}
//_____________________________________________________________________________
//
bool XY::CalculateField(int idx)
{
   if (!X::CalculateField(idx)) return false;
   if (fdC2p[idx]==0 || fdC2m[idx]==0) return false;

   if (idx%(fN1*fN2)<fN1) // C2 lower boundary
      fE2[idx]=(fV[idx]-fV[idx+fN1])/fdC2p[idx];
   else if (idx%(fN1*fN2)>=fN-fN1) // C2 upper boundary
      fE2[idx]=(fV[idx-fN1]-fV[idx])/fdC2m[idx];
   else { // bulk
      fE2[idx]=(fV[idx-fN1]-fV[idx+fN1])/(fdC2m[idx]+fdC2p[idx]);
   }
   return true;
}
#include <vector>
//_____________________________________________________________________________
//
TGraph* XY::GetFieldLineFrom(double x, double y, bool positive)
{
   const char *name = Form("g%.0f%.0f%d",100+x/mm,100+y/mm,positive);
   TGraph *gl=0;
   if (fEgraphs->GetListOfGraphs()) {
      gl = (TGraph*) (fEgraphs->GetListOfGraphs()->FindObject(name));
      if (gl) return gl; // return old graph if it exists
   }

   gl = new TGraph; gl->SetName(name); // create a new graph with a unique name

   int i=0;
   while (true) { // if (x,y) is in crystal
      gl->SetPoint(i,x/cm,y/cm); // add a point to the graph
      if (x>=fC1[fN-1]||x<=fC1[0]||y>=fC2[fN-1]||y<=fC2[0]) {//out of crystal
         if (i==0) // initial point is not in crystal
            Warning("GetFieldLineFrom", "Start point (%.1fcm,%.1fcm)"
                  " is not in crystal! Stop propagating.", x/cm, y/cm);
         break;
      }

      double ex = GetE1(x,y), ey = GetE2(x,y);
      double et = TMath::Sqrt(ex*ex+ey*ey); // total E
      if (et==0) {
         Warning("GetFieldLineFrom", "E@(x=%.1fmm,y=%.1fmm)=%.1fV/cm!",
               x/mm, y/mm, et/volt*cm);
         break;
      }
      double weight=5/(et/volt*mm); // propagate more in weaker field
      double dt=10*ns, mu=50000*cm2/volt/sec; // mu is mobility
      double dx=mu*ex*dt*weight, dy=mu*ey*dt*weight;
      if (gDebug>0)
         Printf("%04d x=%.3fmm, y=%.3fmm, Ex=%.2fV/cm, Ey=%.2fV/cm, "
               "weight=%.3f, dx=%.3fmm, dy=%.3fmm", i, x/mm, y/mm, ex/volt*cm,
               ey/volt*cm, weight, dx/mm, dy/mm);
      if (i>2000) {
         Info("GetFieldLineFrom", "Propagated more than 2000 steps. Stop");
         break;
      }
      if (positive) { x+=dx; y+=dy; } else { x-=dx; y-=dy; }
      i++;
   }

   return gl;
}
