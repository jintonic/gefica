#include <stdlib.h>

#include "PointContactRZ.h"
#include "iostream"
#include "Units.h"
using namespace GeFiCa;
using namespace std;
//for the case which point contact may not fall right on the grid, it will case a great different potential on space close to point contact
//the current design does not include when point contact is close to boundary of the whole detector.
void PointContactRZ::BounardaryOnPointcontact()
{
   int index=FindIdx(PointR,PointDepth,0,n2-1);
   cout<<index<<" "<<fC1[index]<<" "<<fC2[index]<<endl;
   int idxZ=index/n1*n1-1;
   int idxPos=index%n1-1;
   cout<<idxPos<<" "<<fC1[idxPos]<<" "<<fC2[idxPos]<<endl;
   cout<<idxZ<<" "<<fC1[idxZ]<<" "<<fC2[idxZ]<<endl;

   for (int i=0;i<n1;i++) {
      if (PointDepth==0) break;
      fC2[idxZ+i]=PointDepth;//set z for the line where z is closest to pc's depth
      //cout<<idxZ+i<<" "<< fC2[idxZ+i]<<"| ";

      //set steps from depthline
      fDistanceToLeft[idxZ+i]=fC2[idxZ+i]-fC2[idxZ+i-n1];
      fDistanceToRight[idxZ+i]=fC2[idxZ+i+n1]-fC2[idxZ+i];

      //set steps for two line aside depth of pc
      fDistanceToLeft[idxZ+i+n1]=fDistanceToRight[idxZ+i];
      fDistanceToRight[idxZ+i-n1]=fDistanceToLeft[idxZ+i];
   }

   for(int i=0;i<n2;i++) {
      fC1[idxPos+i*n1]=PointR;//set r for the line where r is closest to pc's R 

      //set steps for Radius line
      fDistanceToPrevious[idxPos+i*n1]=fC1[idxPos+i*n1]-fC1[idxPos+i*n1-1];
      fDistanceToNext[idxPos+i*n1]=fC1[idxPos+i*n1+1]-fC1[idxPos+i*n1];

      //set steps for two line aside the previous line
      fDistanceToPrevious[idxPos+i*n1+1]=fDistanceToNext[idxPos+i*n1];
      fDistanceToNext[idxPos+i*n1-1]=fDistanceToPrevious[idxPos+i*n1];
   }
   int idxNeg=n1-idxPos-1;
   cout<<idxNeg<<" "<<fC1[idxNeg]<<" "<<fC2[idxNeg]<<endl;
   for(int i=0;i<n2;i++) {
      fC1[idxNeg+i*n1]=-PointR;//set r for the line where r is closest to pc's R 

      //set steps for Radius line
      fDistanceToPrevious[idxNeg+i*n1]=fC1[idxNeg+i*n1]-fC1[idxNeg+i*n1-1];
      fDistanceToNext[idxNeg+i*n1]=fC1[idxNeg+i*n1+1]-fC1[idxNeg+i*n1];

      //set steps for two line aside the previous line
      fDistanceToPrevious[idxNeg+i*n1+1]=fDistanceToNext[idxNeg+i*n1];
      fDistanceToNext[idxNeg+i*n1-1]=fDistanceToPrevious[idxNeg+i*n1];
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
   BounardaryOnPointcontact();

   // set initial potential values
   for(int i=n;i-->0;) {
      fPotential[i]=(V0+V1)/2;
      // set potential for inner electrodes
      if(fC1[i]>=PointBegin-steplength1/2&&fC1[i]<=PointEnd+steplength1/2&&fC2[i]<=PointDepth+steplength2/2) {
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


}
//_____________________________________________________________________________
//
bool PointContactRZ::CalculateField(EMethod method)
{
   if (!fIsLoaded) Initialize();
   return RZ::CalculateField(method);
}
