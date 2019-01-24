#include <stdlib.h>

#include "PointContactRZ.h"
#include "iostream"
#include "Units.h"
using namespace GeFiCa;
using namespace std;
//for the case which point contact may not fall right on the grid, it will case
//a great different potential on space close to point contact the current
//design does not include when point contact is close to boundary of the whole
//detector.
void PointContactRZ::BoundaryOnPointcontact()
{
   int index=FindIdx(PointR,PointDepth,0,n2-1);
   //cout<<index<<" "<<fC1[index]<<" "<<fC2[index]<<endl;
   int idxZ=(index/n1-1)*n1;
   int idxPos=index%n1-1;
   //cout<<idxPos<<" "<<fC1[idxPos]<<" "<<fC2[idxPos]<<endl;
   //cout<<idxZ<<" "<<fC1[idxZ]<<" "<<fC2[idxZ]<<endl;

   for (int i=0;i<n1;i++) {
      if (PointDepth==0) break;
      fC2[idxZ+i]=PointDepth;//set z for the line where z is closest to pc's depth
      //cout<<idxZ+i<<" "<< fC2[idxZ+i]<<"| ";

      //set steps from depthline
      fdC2m[idxZ+i]=fC2[idxZ+i]-fC2[idxZ+i-n1];
      fdC2p[idxZ+i]=fC2[idxZ+i+n1]-fC2[idxZ+i];

      //set steps for two line aside depth of pc
      fdC2m[idxZ+i+n1]=fdC2p[idxZ+i];
      fdC2p[idxZ+i-n1]=fdC2m[idxZ+i];
   }

   for(int i=0;i<n2;i++) {
      fC1[idxPos+i*n1]=PointR;//set r for the line where r is closest to pc's R 

      //set steps for Radius line
      fdC1m[idxPos+i*n1]=fC1[idxPos+i*n1]-fC1[idxPos+i*n1-1];
      fdC1p[idxPos+i*n1]=fC1[idxPos+i*n1+1]-fC1[idxPos+i*n1];

      //set steps for two line aside the previous line
      fdC1m[idxPos+i*n1+1]=fdC1p[idxPos+i*n1];
      fdC1p[idxPos+i*n1-1]=fdC1m[idxPos+i*n1];
   }
   int idxNeg=n1-idxPos-1;
   //cout<<idxNeg<<" "<<fC1[idxNeg]<<" "<<fC2[idxNeg]<<endl;
   for(int i=0;i<n2;i++) {
      fC1[idxNeg+i*n1]=-PointR;//set r for the line where r is closest to pc's R 

      //set steps for Radius line
      fdC1m[idxNeg+i*n1]=fC1[idxNeg+i*n1]-fC1[idxNeg+i*n1-1];
      fdC1p[idxNeg+i*n1]=fC1[idxNeg+i*n1+1]-fC1[idxNeg+i*n1];

      //set steps for two line aside the previous line
      fdC1m[idxNeg+i*n1+1]=fdC1p[idxNeg+i*n1];
      fdC1p[idxNeg+i*n1-1]=fdC1m[idxNeg+i*n1];
   }
}
void PointContactRZ::BoundaryonWarpAround()
{
   int index=FindIdx(ContactInnerR,0,0,n2-1);
   if (index>n1)index-=n1;
   for(int i=index;i<n;i+=n1)
   {
      fC1[i]=ContactInnerR;
      fdC1m[i]=fC1[i]-fC1[i-1];
      fdC1p[i]=fC1[i+1]-fC1[i];
      fdC1m[i+1]=fdC1p[i];
      fdC1p[i-1]=fdC1m[i];

   }

   index=FindIdx(-ContactInnerR,0,0,n2-1)-1;
   if (index>n1)index-=n1;

   for(int i=index;i<n;i+=n1)
   {
      fC1[i]=-ContactInnerR;
      fdC1m[i]=fC1[i]-fC1[i-1];
      fdC1p[i]=fC1[i+1]-fC1[i];
      fdC1m[i+1]=fdC1p[i];
      fdC1p[i-1]=fdC1m[i];

   }
}
void PointContactRZ::Initialize()
{
   if (n1%2==1) {
         Error("Initialize", "Number of grids in R cannot be odd, abort!");
         abort();
   }

   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (ZLowerBound>=ZUpperBound) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            ZLowerBound, ZUpperBound);
      return;
   }
   double RUpperBound,RLowerBound,PointBegin,PointEnd;
   RUpperBound=Radius;
   RLowerBound=-Radius;
   PointBegin=-PointR;
   PointEnd=PointR;
   double steplength1=(RUpperBound-RLowerBound)/(n1-1);
   double steplength2=(ZUpperBound-ZLowerBound)/(n2-1);
   SetStepLength(steplength1,steplength2);
   for(int i=n;i-->0;) 
   {
      fC1[i]=fC1[i]+RLowerBound;
   } 
   BoundaryOnPointcontact();
   BoundaryonWarpAround();

   // set initial potential values
   for(int i=n;i-->0;) {
      fPotential[i]=(V0+V1)/2;//common this line for finding depleat voltage
      // set potential for inner electrodes
      if(fC1[i]>=PointBegin&&fC1[i]<=PointEnd&&fC2[i]<=PointDepth) {
         fPotential[i]=V1;
         fIsFixed[i]=true;
      }
   }
   // set potential for outer electrodes
   for(int i=n-1;i>=n-n1;i--) {
      fIsFixed[i]=true;
      fPotential[i]=V0;
   }
   for(int i=0;i<n-n1;i=i+n1) {
      fIsFixed[i]=true;
      fIsFixed[i+n1-1]=true;
      fPotential[i]=V0;
      fPotential[i+n1-1]=V0;
   }
   for (int i=0;i<n1;i++)
   {
      if(fC1[i]>=ContactInnerR||fC1[i]<=-ContactInnerR)
      {
         fIsFixed[i]=true;
         fPotential[i]=V0;
      }
   }
}
//_____________________________________________________________________________
//
bool PointContactRZ::CalculatePotential(EMethod method)
{
   if (!fIsLoaded) Initialize();
// this commentd block are slow depletion voltage finder
 while(0)
 {
    RZ::CalculatePotential(method);
    if(!X::IsDepleted())
    {
        int maxn=Findmax();
        int minn=Findmin();
        if(V0>V1)
        {
           V0=(fPotential[maxn]-V0)*1.01+V0;
           V1=-(V1-fPotential[minn])*1.01+V1;
        }
        else
        {
           V1=(fPotential[maxn]-V1)*1.01+V1;
           V0=-(V0-fPotential[minn])*1.01+V0;
        }
        Initialize();
        RZ::CalculatePotential(method);
        cout<<V0<<" "<<V1<<endl;
    }
    else break;
 }
   return RZ::CalculatePotential(method);
}
//_____________________________________________________________________________
//
bool PointContactRZ::CalculateField(int idx)
{
   if (!XY::CalculateField(idx)) return false;

   if (fC2[idx]>PointDepth-fdC2m[idx]
         && fC2[idx]<PointDepth+fdC2p[idx]) // PC top boundary
      fE2[idx]=(fPotential[idx]-fPotential[idx+n1])/fdC2p[idx];
   if (fC1[idx]>-PointR-fdC1m[idx]
         && fC1[idx]<-PointR+fdC1p[idx]) // PC left boundary
      fE1[idx]=(fPotential[idx]-fPotential[idx-1])/fdC1m[idx];
   if (fC1[idx]>PointR-fdC1m[idx]
         && fC1[idx]<PointR+fdC1p[idx]) // PC right boundary
      fE1[idx]=(fPotential[idx]-fPotential[idx+1])/fdC1p[idx];

   return true;
}
#include <fstream>
bool PointContactRZ::SaveFieldasFieldgen(const char * fout)
{
   ofstream outfile(fout);

   outfile<<"# height "<< ZUpperBound-ZLowerBound;        
   outfile<<"\n# xtal_radius "<<Radius;
   outfile<<"\n# pc_length   "<<PointDepth;        
   outfile<<"\n# pc_radius   "<<PointR;         
   outfile<<"\n# wrap_around_radius "<<ContactInnerR; 
   outfile<<"\n# grid size on r "<<fdC1p[0];
   outfile<<"\n# grid size on z "<<fdC2p[0];
   outfile<<"\n# impurity_z0  "<<fImpurity[0];
   outfile<<"\n# xtal_HV      "<<V1;
   outfile<<"\n# max_iterations "<<MaxIterations;
   outfile<<"\n# ";
   outfile<<"\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)";
   for (int i=0;i<n;i++) {
      double E=sqrt(fE1[i]*fE1[i]+fE2[i]*fE2[i]);
      outfile<<"\n"<<fC1[i]<<"  "<<fC2[i]<<"  "<<fPotential[i]<<"  "<<E<<"  "<<fE1[i]<<"  "<<fE2[i];
   }
   outfile.close();
   return true;
}
