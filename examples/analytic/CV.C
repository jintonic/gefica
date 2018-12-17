// create C-V curve for an ideal planar detector
#include<math.h>
// definition of necessary units
static const double cm=1;
static const double cm3=cm*cm*cm;
static const double volt=1;
static const double C=1; // Coulomb
static const double e=1.6e-19*C; // elementary charge
static const double epsilon0=8.854187817e-14*C/volt/cm; // vacuum permittivity
// https://link.springer.com/chapter/10.1007/10832182_519
static const double epsilonGe=15.8; // Ge dielectric constant
double VtoC(double V,double totalD)
{
   double impurity=-1e10*e/cm3;
   double d=sqrt(-2*epsilonGe*epsilon0*V/impurity);
   //V=ax^2+c2x+c1
   //c1=0, when V=0 at x=0
   //c2=-ad, when E=dV/dx=0 at x=0, where just depleted
   //a is impurity/epsilon
   //V=-ax^2/2, solve V when x=d
   if (d>totalD)d=totalD;
   //when it depleted whole detector
   double Cap=epsilon0*epsilonGe/d;
   return Cap;

}
void CV()
{
   int n=200;
   double Cap[n],V[n];
   for (int i=1;i<=n;i++)
   {
      V[i-1]=i*10*volt;
      Cap[i-1]=VtoC(V[i-1],1*cm);
   }
   TGraph *gr  = new TGraph(n,V,Cap);
   TCanvas *c1 = new TCanvas("c1","Graph Draw Options",
                                         200,10,600,400);

   // draw the graph with axis, continuous line, and put
   // a * at each point
   gr->Draw("AC*");

}
