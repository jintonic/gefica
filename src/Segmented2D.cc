#include <stdlib.h>

#include "Segmented2D.h"
#include "Units.h"
using namespace GeFiCa;
using namespace std;
void Segmented2D::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (RLowerBound>=RUpperBound) {
      Warning("Initialize",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            RLowerBound, RUpperBound);
      return;
   }
   if(SegmentNum==0)
   {
      Warning("Initialize",
            "total Segment number cannot be zero!");
      return;
   }
   if (SegmentNum<0)SegmentNum=-SegmentNum;
   if (SegmentID<0)SegmentID=-SegmentID;
   if(SegmentID>SegmentNum)SegmentID=SegmentID%SegmentNum;

   double steplength1=(RUpperBound-RLowerBound)/(n1-1);
   double steplength2=2*TMath::Pi()/(n2);
   double SegmentUpperBound=2*TMath::Pi()*SegmentID/SegmentNum;
   double SegmentLowerBound=2*TMath::Pi()*(SegmentID-1)/SegmentNum;
   SetStepLength(steplength1,steplength2);
   for(int i=0;i<n;i++)
   {
      fC1[i]+=RLowerBound;
      if(i%n1==0)
      {
         fIsFixed[i]=true;
         if(SegmentID==0)fPotential[i]=V1;
         else fPotential[i]=V0;
      }
      else if(i%n1==n1-1)
         //need shift boundary as pointcontact
      {
         fIsFixed[i]=true;
         if(SegmentID==0)fPotential[i]=V0;
         else if(fC2[i]<=SegmentUpperBound&&fC2[i]>=SegmentLowerBound)
         {
            fPotential[i]=V1;
         }
         else if(fC2[i]>=SegmentUpperBound||fC2[i]<=SegmentLowerBound)
         {
            fPotential[i]=V0;
         }
      }
      else 
      {
         fIsFixed[i]=false;
         fPotential[i]=0;
      }


   }
}
//_____________________________________________________________________________
//
bool Segmented2D::CalculatePotential(EMethod method)
{
   if (!fIsLoaded) Initialize();
   return RhoPhi::CalculatePotential(method);
}
