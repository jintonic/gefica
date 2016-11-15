#include "Planar1D.h"
using namespace GeFiCa;

void Planar1D::Initialize()
{
   // The step length is calculated with the following equation:
   // BEGIN_HTML
   // <pre>
   //      double stepLength=(UpperBound-LowerBound)/(n-1);
   // </pre>
   // END_HTML
   // If the inner radius is not larger than the outer radius,
   // no grid will be created
   if (LowerBound>=UpperBound) {
      Warning("CreateGridWithFixedStepLength",
            "Lower bound (%f) >= upper bound (%f)! No grid is created!",
            LowerBound, UpperBound);
      return;
   }
   double steplength=(UpperBound-LowerBound)/(n-1);
   SetStepLength(steplength);
   for (int i=0; i<n; i++) fC1[i]=fC1[i]+LowerBound;
   fIsFixed[0]=true;
   fIsFixed[n-1]=true;
   double slope = (Vpos-Vneg)/(n-1);
   for (int i=0; i<n; i++) fPotential[i]=Vneg+slope*i;
}
//_____________________________________________________________________________
//
//Here we calculate the field analytically using the fomrula (d^2 phi/d phi^2) = rho/epsilon.
//Rho in this case is just the impurity density. This gives us a differential equation 
//the type phi(x)= ax^2 + bx + c . We can find a, b, and c by evaluating this equation
//at the boundary conditions where we know the potential
// we know:
// a=rho/(2*epsilon)
// b=(Vneg-Vpos-a(LowerBound^2-UpperBound^2))/(LowerBound-UpperBound)
// c=Vneg-a*LowerBound^2-LowerBound(Vneg-Vpos-a*(LowerBound^2-UpperBound^2))/(LowerBound-UpperBound)



// bool Planar1D::Analytic()
// {
//    double d=LowerBound-UpperBound;
//    double dSquare=LowerBound*lowerBound-UpperBound*UpperBound;
//    double Aconst=fImpurity[i]*Qe/2/epsilon;
//    double Bconst=(Vneg-Vpos-Aconst*dSquare)/d;
//    double Cconst=Vneg-Aconst*LowerBound*LowerBound-LowerBound*(Vneg-Vpos-Aconst*dSquare)/d;
//    for (int i=0; i<n; i++) {
//      fPotential[i] = Aconst*fC1[i]*fC1[i]+Bconst*fC1[i]+Cconst;
//      fE1[i]=(fPotential[i+1]-fPotential[i-1])/(fDistanceToNext[i]+fDistanceToPrevious[i]);
//    }
//   return true;
// }



bool Planar1D::Analytic()
{
   double d=UpperBound-LowerBound;//thickness or depth of the detector
   double cnst1=fPotential[0];
   double cnst2=(fPotential[n-1]-fImpurity[n-1]*Qe/2/epsilon*d*d-cnst1)/d;
   for (int i=0; i<n; i++) {
      fPotential[i] = fImpurity[i]*Qe/2/epsilon*fC1[i]*fC1[i]+cnst2*fC1[i]+cnst1;
      fE1[i]=(fPotential[i+1]-fPotential[i-1])/(fDistanceToNext[i]+fDistanceToPrevious[i]);
   }
   return true;
}
//_____________________________________________________________________________
//
bool Planar1D::CalculateField(EMethod method)
{
   if(!fIsLoaded)Initialize();
   return X::CalculateField(method);
}
