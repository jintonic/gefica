#include <cmath>
#include <vector>
using namespace std;

#include "XY.h"
#include "Units.h"
using namespace GeFiCa;

void XY::OverRelaxAt(size_t idx)
{
   if (fIsFixed[idx])return;

   double vym,vyp,vxm,vxp;
   if (idx>=N1) vym=V[idx-N1];
   else vym=V[idx];
   if (idx>=V.size()-N1) vyp=V[idx];
   else vyp=V[idx+N1];
   if (idx%N1==0) vxm=V[idx];
   else vxm=V[idx-1];
   if (idx%N1==N1-1) vxp=V[idx];
   else vxp=V[idx+1];

   double tmp = (Src[idx]
         + (vxp/dC1m[idx]+vxm/dC1p[idx])*2/(dC1m[idx]+dC1p[idx])
         + (vyp/dC2p[idx]+vym/dC2m[idx])*2/(dC2p[idx]+dC2m[idx]))/
      ((1/dC1m[idx]+1/dC1p[idx])*2/(dC1m[idx]+dC1p[idx])
       + (1/dC2p[idx]+1/dC2m[idx])*2/(dC2p[idx]+dC2m[idx]));
   
   tmp=RelaxationFactor*(tmp-V[idx])+V[idx];

   double min=vxm, max=vxm;
   if (min>vxp) min=vxp;
   if (min>vyp) min=vyp;
   if (min>vym) min=vym;
   
   if (max<vxp) max=vxp;
   if (max<vyp) max=vyp;
   if (max<vym) max=vym;
   
   if (tmp<min) {
      V[idx]=min;
      fIsDepleted[idx]=false;
   } else if (tmp>max) {
      V[idx]=max;
      fIsDepleted[idx]=false;
   } else
      fIsDepleted[idx]=true;

   if(fIsDepleted[idx]||V.begin()==V.end()) V[idx]=tmp;
}
//_____________________________________________________________________________
//
double XY::GetData(const vector<double>& data,
      double x, double y, double z) const
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
   size_t idx=GetIdxOfPointNear(x,y); // always exists

   
   bool tl=false; // existence of top left grid point
   bool br=false; // existence of bottom right grid point
   bool bl=false; // existence of bottom left grid point
   if (idx%N1!=0) tl=true; // not left boundary
   if (idx>=N1) br=true; // not bottom boundary
   if (tl&&bl) bl=true; // neither left nor bottom boundary
//
   if(!tl&&!br&&!bl) {
      return data[idx];
   } else if(tl&&!br&&!bl) {
      return twopoint({data[idx-1],data[idx]},{x,y},{C1[idx-1],C1[idx]});
   }
   else if(!tl&&!bl&&br)
   {
      return twopoint({data[idx-N1],data[idx]},{x,y},{C2[idx-N1],C2[idx]});
   }
   else
   {
      //no bound case
      if(dC1p[idx-1]==dC1m[idx]&&dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&dC2p[idx-N1]==dC2m[idx])
      {
         return fourpoint({data[idx-1],data[idx],data[idx-N1-1],data[idx-N1]},{x,y},{C1[idx-1],C1[idx],C1[idx-N1-1],C1[idx-N1]},{C2[idx-1],C2[idx],C2[idx-N1-1],C2[idx-N1]});
      }
      //topleft case o
      if(dC1p[idx-1]!=dC1m[idx]&&
            dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]!=dC2m[idx-1]&&
            dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint({tmv,data[idx],bmv,data[idx-N1]},{x,y},{xb,C1[idx],xb,C1[idx-N1]},{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint({mlv,mmv,data[idx-N1-1],bmv},{x,y},{C1[idx-1],xb,C1[idx-1],xb},{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint({tmv,mmv,mlv},{x,y},{xb,xb,C1[idx-1]},{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint({tmv,tmv,mlv},{x,y},{xb,C1[idx-1],C1[idx-1]},{C2[idx],C2[idx],yb});
         }
      }
      //topright case o
      if(dC1p[idx-1]!=dC1m[idx]&&
            dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&
            dC2p[idx-N1]!=dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint({tmv,data[idx],bmv,data[idx-N1]},{x,y},{xb,C1[idx],xb,C1[idx-N1]},{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint({mlv,mmv,data[idx-N1-1],bmv},{x,y},{C1[idx-1],xb,C1[idx-1],xb},{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint({tmv,mmv,mlv},{x,y},{xb,xb,C1[idx-1]},{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint({tmv,tmv,mlv},{x,y},{xb,C1[idx-1],C1[idx-1]},{C2[idx],C2[idx],yb});
         }
      }
      if(dC1p[idx-1]==dC1m[idx]&&dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint({tmv,data[idx],bmv,data[idx-N1]},{x,y},{xb,C1[idx],xb,C1[idx-N1]},{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint({mlv,mmv,data[idx-N1-1],bmv},{x,y},{C1[idx-1],xb,C1[idx-1],xb},{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint({tmv,mmv,mlv},{x,y},{xb,xb,C1[idx-1]},{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint({tmv,tmv,mlv},{x,y},{xb,C1[idx-1],C1[idx-1]},{C2[idx],C2[idx],yb});
         }
      }
      if(dC1p[idx-1]==dC1m[idx]&&dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint({tmv,data[idx],bmv,data[idx-N1]},{x,y},{xb,C1[idx],xb,C1[idx-N1]},{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint({mlv,mmv,data[idx-N1-1],bmv},{x,y},{C1[idx-1],xb,C1[idx-1],xb},{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint({tmv,mmv,mlv},{x,y},{xb,xb,C1[idx-1]},{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint({tmv,tmv,mlv},{x,y},{xb,C1[idx-1],C1[idx-1]},{C2[idx],C2[idx],yb});
         }
      }
      if(dC1p[idx-1]==dC1m[idx]&&dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint({tmv,data[idx],bmv,data[idx-N1]},{x,y},{xb,C1[idx],xb,C1[idx-N1]},{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint({mlv,mmv,data[idx-N1-1],bmv},{x,y},{C1[idx-1],xb,C1[idx-1],xb},{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint({tmv,mmv,mlv},{x,y},{xb,xb,C1[idx-1]},{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint({tmv,tmv,mlv},{x,y},{xb,C1[idx-1],C1[idx-1]},{C2[idx],C2[idx],yb});
         }
      }
      if(dC1p[idx-1]==dC1m[idx]&&dC1p[idx-N1-1]==dC1m[idx-N1]&&
            dC2p[idx-N1-1]==dC2m[idx-1]&&dC2p[idx-N1]==dC2m[idx])
      {
         double xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
         double yb=dC2m[idx-1]<dC2p[idx-N1-1] ? C2[idx-1]-dC2m[idx-1] : C2[idx-N1-1]+dC2p[idx-N1-1];
         double bmv=(C1[idx-N1]-xb)/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1-1]+(xb-C1[idx-N1-1])/(C1[idx-N1]-C1[idx-N1-1])*data[idx-N1];
         double tmv=fIsFixed[idx]?data[idx]:data[idx-1];
         double mlv=fIsFixed[idx-1]?data[idx-1]:data[idx-N1-1];
         double mmv=(C2[idx]-yb)/(C2[idx]-C2[idx-N1])*bmv+(yb-C2[idx-N1])/(C2[idx]-C2[idx-N1])*tmv;
         if(x>xb)
         {
            return fourpoint({tmv,data[idx],bmv,data[idx-N1]},{x,y},{xb,C1[idx],xb,C1[idx-N1]},{C2[idx],C2[idx],C2[idx-N1],C2[idx-N1]});
         }
         else if(y<yb)
         {
            return fourpoint({mlv,mmv,data[idx-N1-1],bmv},{x,y},{C1[idx-1],xb,C1[idx-1],xb},{yb,yb,C2[idx-N1],C2[idx-N1]});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb>y))
         {
            return threepoint({tmv,mmv,mlv},{x,y},{xb,xb,C1[idx-1]},{C2[idx],yb,yb});
         }
         else if((C2[idx]-yb)/(xb-C1[idx-1]*(x-C1[idx-1])+yb<y))
         {
            return threepoint({tmv,tmv,mlv},{x,y},{xb,C1[idx-1],C1[idx-1]},{C2[idx],C2[idx],yb});
         }
      }
   }
//   double tmv; // interpolated value at (x, C2[idx])
//   double bmv; // interpolated value at (x, C2[idx-N1])
//   double dyp; // distance between (x,y) and (x, C2[idx]) or a boundary point
//   double dym; // distance between (x,y) and (x, C2[idx-N1]) or a boundary point
//
//   if (tl) { // interpolate tl & tr points if both tl and tr exist
//      double trv,tlv,aa,ab;
//      double xb; // x of boundary crossing point on horizontal lines
//      // in case of normal condistion xb == C1[idx]
//      xb=dC1m[idx]<dC1p[idx-1] ? C1[idx]-dC1m[idx] : C1[idx-1]+dC1p[idx-1];
//      if (x<xb) { // xb is on the right of x
//         tlv=data[idx-1];
//         if(fIsFixed[idx-1]) trv=data[idx-1]; // left is outside of crystal
//         else trv=data[idx]; // right is outside of crystal or normal condition
//         aa=x-C1[idx-1];
//         ab=xb-x;
//      } else { // xb is on the left of x
//         trv=data[idx];
//         if (fIsFixed[idx]) tlv=data[idx]; // right is outside of crystal
//         else tlv=data[idx-1]; // left is outside of crystal or normal 
//         ab=C1[idx]-x;
//         aa=x-xb;
//      }
//      tmv = trv*aa/(aa+ab) + tlv*ab/(aa+ab);
//   } else tmv=data[idx];
//
//   if(br&&bl) {
//      double brv,blv,aa,ab;
//      double xb; // x of boundary crossing point on horizontal lines
//      // in case of normal condistion xb == C1[idx]
//      xb=dC1m[idx-N1]<dC1p[idx-N1-1] ? C1[idx-N1]-dC1m[idx-N1] : C1[idx-N1-1]+dC1p[idx-N1-1];
//      if (x<xb) { // xb is on the right of x
//         blv=data[idx-N1-1];
//         if(fIsFixed[idx-N1-1]) brv=data[idx-N1-1]; // left is outside of crystal
//         else brv=data[idx-N1]; // right is outside of crystal or normal condition
//         aa=x-C1[idx-N1-1];
//         ab=xb-x;
//      } else { // xb is on the left of x
//         brv=data[idx-N1];
//         if (fIsFixed[idx-N1]) blv=data[idx-N1]; // right is outside of crystal
//         else blv=data[idx-N1-1]; // left is outside of crystal or normal 
//         ab=C1[idx-N1]-x;
//         aa=x-xb;
//      }
//      bmv = brv*aa/(aa+ab) + blv*ab/(aa+ab);
//   } else if (br==false && bl==false) {
//      bmv = tmv;
//   } else
//      bmv = data[idx-N1];
//        
//   double ab=(-x+C1[idx])/dC1m[idx];
//   double aa=1-ab;
//   double ba=(-y+C2[idx])/dC2m[idx];
//   double bb=1-ba;
//
//   if (gDebug>0) {
//      Info("GetData","C1[%d]=%f, ",idx,C1[idx]);
//   }
//   return (tlv*ab+trv*aa)*bb+(blv*ab+brv*aa)*ba;
}
//_____________________________________________________________________________
//
bool XY::CalculateField(size_t idx)
{
   //if (!X::CalculateField(idx)) return false;
   if (dC2p[idx]==0 || dC2m[idx]==0) return false;

   if (idx%(N1*N2)<N1) // C2 lower boundary
      E2[idx]=(V[idx]-V[idx+N1])/dC2p[idx];
   else if (idx%(N1*N2)>=V.size()-N1) // C2 upper boundary
      E2[idx]=(V[idx-N1]-V[idx])/dC2m[idx];
   else { // bulk
      E2[idx]=(V[idx-N1]-V[idx+N1])/(dC2m[idx]+dC2p[idx]);
   }
   return true;
}
//_____________________________________________________________________________
//
FieldLine* XY::GetFieldLineFrom(double x, double y, bool positive)
{
   FieldLine *fl = new FieldLine;

   size_t i=0;
   while (true) { // if (x,y) is in crystal
      fl->C1.push_back(x/cm); fl->C2.push_back(y/cm);
      if (x>=C1[V.size()-1]||x<=C1[0]||y>=C2[V.size()-1]||y<=C2[0]) {//out of crystal
         if (i==0) // initial point is not in crystal
            Warning("GetFieldLineFrom", "Start point (%.1fcm,%.1fcm)"
                  " is not in crystal! Stop propagating.", x/cm, y/cm);
         break;
      }

      double ex = GetE1(x,y), ey = GetE2(x,y);
      double et = sqrt(ex*ex+ey*ey); // total E
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

   return fl;
}

double XY::twopoint(dataset[2],tarlocationset[1],pointxset[2])
{
   double ab=abs(pointxset[0]-tarlocationset[0])/abs(pointxset[1]-pointxset[0]);
   double aa=1-ab;
   return dataset[0]*aa/(aa+ab) + dataset[1]*ab/(aa+ab);
}

double XY::threepoint(dataset[3],tarlocationset[2],pointxset[3],pointyset[3])
{
   double x=tarlocationset[0];
   double x1=pointxset[0];
   double x2=pointxset[1];

   double y=tarlocationset[1];
   double y1=pointyset[0];
   double y2=pointyset[1];
   double y3=pointyset[2];
   //x=(1-u-v)*x1+u*x2+v*x3
   //y=(1-u-v)*y1+u*y2+v*y3
   double v=((y2-y1)*(x-x1)-(y-y1)*(x2-x1))/((x3-x1*(y2-y1)-(y3-y1)(x2-x1)));
   double u=(y-y1-v*(y3-y1))/(y2-y1);

   return dataset[0]*(1-u-v)+dataset[1]*u+dataset[2]*v;
}
//0---------1
//|         |
//|         |
//|         |
//|         |
//|         |
//2---------3
double XY::fourpoint(dataset[4],tarlocationset[2],pointxset[4],pointyset[4])
{
   double ab=(pointxset[1]-tarlocationset[0])/(pointxset[1]-pointxset[0]);
   double aa=1-ab;
   double ba=(tarlocationset[1]-pointyset[1])/(pointyset[3]-pointyset[1]);
   double bb=1-ba;
   
   return (dataset[1]*ab+dataset[0]*aa)*bb+(dataset[2]*ab+dataset[3]*aa)*ba;
}
