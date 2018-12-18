#include <stdlib.h>

#include "Siegfried2D.h"
#include "Units.h"
using namespace GeFiCa;
using namespace std;
void Siegfried2D::Initialize()
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
   double steplength1=(RUpperBound-RLowerBound)/(n1-1);
   double steplength2=2*TMath::Pi()/(n2);
   SetStepLength(steplength1,steplength2);
   for(int i=0;i<n;i++)
   {
      fC1[i]+=RLowerBound;
      if(i%n1==0)
      {
         fIsFixed[i]=true;
         fPotential[i]=V0;
      }
      else if(i%n1==n1-1&&fC2[i]<=SegmentSize)
         //need shift boundary as pointcontact
      {
         fIsFixed[i]=true;
         fPotential[i]=V1;
      }
      else if(i%n1==n1-1&&fC2[i]>SegmentSize)
      {
         fIsFixed[i]=true;
         fPotential[i]=V0;
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
bool Siegfried2D::CalculatePotential(EMethod method)
{
   if (!fIsLoaded) Initialize();
   return RhoPhi::CalculatePotential(method);
}
