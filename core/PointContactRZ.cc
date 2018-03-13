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
    int StartinR=index/n1*n1;
    int StartinZ=index%n1;
    cout<<StartinZ<<" "<<fC1[StartinZ]<<" "<<fC2[StartinZ]<<endl;
    cout<<StartinR<<" "<<fC1[StartinR]<<" "<<fC2[StartinR]<<endl;

    for (int i=0;i<n1;i++)
    {
        fC2[StartinR+i]=PointDepth;//set z for the line where z is closest to pc's depth
        //cout<<StartinR+i<<" "<< fC2[StartinR+i]<<"| ";

        //set steps from depthline
        fDistanceToPrevious[StartinR+i]=fC2[StartinR+i]-fC2[StartinR+i-n1];
        fDistanceToNext[StartinR+i]=fC2[StartinR+i+n1]-fC2[StartinR+i];

        //set steps for two line aside depth of pc
        fDistanceToPrevious[StartinR+i+n1]=fDistanceToNext[StartinR+i];
        fDistanceToNext[StartinR+i-n1]=fDistanceToPrevious[StartinR+i];


    }

    for(int i =0;i<n2;i++)
    {
        fC1[StartinZ+i*n1]=PointR;//set r for the line where r is closest to pc's R 

        //set steps for Radius line
        fDistanceToLeft[StartinZ+i*n1]=fC1[StartinZ+i*n1]-fC1[StartinZ+i*n1-1];
        fDistanceToRight[StartinZ+i*n1]=fC1[StartinZ+i*n1+1]-fC1[StartinZ+i*n1];

        //set steps for two line aside the previous line
        fDistanceToLeft[StartinZ+i*n1+1]=fDistanceToRight[StartinZ+i*n1];
        fDistanceToRight[StartinZ+i*n1-1]=fDistanceToLeft[StartinZ+i*n1];
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
