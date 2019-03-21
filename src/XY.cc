#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "XY.h"
#include "Units.h"
using namespace GeFiCa;

XY::XY(int nx, int ny, const char *name, const char *title)
   : X(nx*ny, name, title)
{
   N1=nx; // N1 is set to nx*ny through X constructor, it is fixed here
   N2=ny;
   fEgraphs = new TMultiGraph("gEs","electric field lines");
}
//_____________________________________________________________________________
//
XY::~XY()
{
   if (E2) delete[] E2;
   if (C2) delete[] C2;
   if (dC2m) delete[] dC2m;
   if (dC2p) delete[] dC2p;
   if (fEgraphs) delete fEgraphs;
}
//_____________________________________________________________________________
//
void XY::SetStepLength(double steplength1,double steplength2)
{
   //set field step length
   X::SetStepLength(steplength1);
   for (int i=0;i<fN;i++) {
      if(i>N1-1)C2[i]=C2[i-N1]+steplength2;
      else C2[i]=0;
      if(i%N1==0)C1[i]=0;
      else C1[i]=C1[i-1]+steplength1;

      E2[i]=0;
      dC2m[i]=steplength2;
      dC2p[i]=steplength2;
   }
}
//_____________________________________________________________________________
//
#include <iostream>
using namespace std;
void XY::OverRelaxAt(int idx)
{
   // 2nd-order Runge-Kutta Successive Over-Relaxation
   if (fIsFixed[idx])return;
   double density=fImpurity[idx]*Qe;
   double h2=dC1m[idx];
   double h3=dC1p[idx];
   double h4=dC2m[idx];
   double h1=dC2p[idx];
   double pym,pyp,pxm,pxp;
   if(idx>=N1)pym=V[idx-N1];
   else pym=V[idx];
   if(idx>=fN-N1)pyp=V[idx];
   else pyp=V[idx+N1];
   if(idx%N1==0)pxm=V[idx];
   else pxm=V[idx-1];
   if(idx%N1==N1-1)pxp=V[idx];
   else pxp=V[idx+1];
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
double XY::GetData(double x, double y, double z, double *data)
{
   // https://codeplea.com/triangular-interpolation
   //       tmv
   // +-aa--+--ab--+(C1[idx], C2[idx])
   // |     ^      |
   // |  dyp|      ba
   // |     |(x,y) |
   // +<----+----->+
   // |     |      |
   // |  dym|      bb
   // |     v      |
   // +-----+------+
   //       bmv
   int idx=FindIdx(x,y); // always exists

   bool tl=false; // existence of top left grid point
   bool br=false; // existence of bottom right grid point
   bool bl=false; // existence of bottom left grid point
   if (idx%N1!=0) tl=true; // not left boundary
   if (idx>=N1) br=true; // not bottom boundary
   if (tl&&bl) bl=true; // neither left nor bottom boundary

   double tmv; // interpolated value at (x, C2[idx])
   double bmv; // interpolated value at (x, C2[idx-N1])
   double dyp; // distance between (x,y) and (x, C2[idx]) or a boundary point
   double dym; // distance between (x,y) and (x, C2[idx-N1]) or a boundary point

   if (tl) { // interpolate tl & tr points if both tl and tr exist
      double trv,tlv,aa,ab;
      double xb; // x of boundary crossing point on horizontal lines
      // in case of normal condistion xb == C1[idx]
      xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
      if (x<xb) { // xb is on the right of x
         tlv=data[idx-1];
         if(fIsFixed[idx-1]) trv=data[idx-1]; // left is outside of crystal
         else trv=data[idx]; // right is outside of crystal or normal condition
         aa=x-C1[idx-1];
         ab=xb-x;
      } else { // xb is on the left of x
         trv=data[idx];
         if (fIsFixed[idx]) tlv=data[idx]; // right is outside of crystal
         else tlv=data[idx-1]; // left is outside of crystal or normal 
         ab=C1[idx]-x;
         aa=x-xb;
      }
      tmv = trv*aa/(aa+ab) + tlv*ab/(aa+ab);
   } else tmv=data[idx];

   if(br&&bl) {
      double brv,blv,aa,ab;
      double xb; // x of boundary crossing point on horizontal lines
      // in case of normal condistion xb == C1[idx]
      xb=dC1m[idx-N1]<dC1p[idx-N1-1] ? C1[idx-N1]-dC1m[idx-N1] : C1[idx-N1-1]+dC1p[idx-N1-1];
      if (x<xb) { // xb is on the right of x
         blv=data[idx-N1-1];
         if(fIsFixed[idx-N1-1]) brv=data[idx-N1-1]; // left is outside of crystal
         else brv=data[idx-N1]; // right is outside of crystal or normal condition
         aa=x-C1[idx-N1-1];
         ab=xb-x;
      } else { // xb is on the left of x
         brv=data[idx-N1];
         if (fIsFixed[idx-N1]) blv=data[idx-N1]; // right is outside of crystal
         else blv=data[idx-N1-1]; // left is outside of crystal or normal 
         ab=C1[idx-N1]-x;
         aa=x-xb;
      }
      bmv = brv*aa/(aa+ab) + blv*ab/(aa+ab);
   } else if (br==false && bl==false) {
      bmv = tmv;
   } else
      bmv = data[idx-N1];
        
   double ab=(-x+C1[idx])/dC1m[idx];
   double aa=1-ab;
   double ba=(-y+C2[idx])/dC2m[idx];
   double bb=1-ba;

   if (gDebug>0) {
      Info("GetData","C1[%d]=%f, ",idx,C1[idx]);
   }
   return (tlv*ab+trv*aa)*bb+(blv*ab+brv*aa)*ba;
}
//_____________________________________________________________________________
//
bool XY::CalculateField(int idx)
{
   if (!X::CalculateField(idx)) return false;
   if (dC2p[idx]==0 || dC2m[idx]==0) return false;

   if (idx%(N1*N2)<N1) // C2 lower boundary
      E2[idx]=(V[idx]-V[idx+N1])/dC2p[idx];
   else if (idx%(N1*N2)>=fN-N1) // C2 upper boundary
      E2[idx]=(V[idx-N1]-V[idx])/dC2m[idx];
   else { // bulk
      E2[idx]=(V[idx-N1]-V[idx+N1])/(dC2m[idx]+dC2p[idx]);
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
      if (x>=C1[fN-1]||x<=C1[0]||y>=C2[fN-1]||y<=C2[0]) {//out of crystal
         if (i==0) // initial point is not in crystal
            Warning("GetFieldLineFrom", "Start point (%.1fcm,%.1fcm)"
                  " is not in crystal! Stop propagating.", x/cm, y/cm);
         break;
      }

      double ex = GetE1(x,y), ey = GetE2(x,y);
      double et = TMath::Sqrt(ex*ex+ey*ey); // total E
      if (et==0) {
         Warning("GetFieldLineFrom", "E@(x=%.1fmm,y=%.1fmm)=%.1V/cm!",
               x/mm, y/mm, et/volt*cm);
         break;
      }
      double weight=5/(et/volt*mm); // propagate more in weaker field
      double dt=10*ns, mu=50000*cm2/volt/sec; // mu is mobility
      double dx=mu*ex*dt*weight, dy=mu*ey*dt*weight;
      if (gDebug>0)
         Printf("%04d x=%.3fmm, y=%.3fmm, Ex=%.2V/cm, Ey=%.2V/cm, "
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
