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
double NumericalFieldGenerate(double V, double totalD)
{
   int n=2001;
   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(n);
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->UpperBound=totalD*GeFiCa::cm;
   detector->V1=0*GeFiCa::volt;
   detector->V0=V*GeFiCa::volt;

   TF1 *im=new TF1("f","-1e10");
   detector->SetImpurity(im);
   bool *connect=new bool[n];
   for(int i=0;i<n;i++)
   {
      connect[i]=false;
   }
   connect[0]=true;

   detector->CalculatePotential(GeFiCa::kSOR2);
   double d=0;
   for(int i=0;i<n-1;i++)
   {
      if(detector->fIsDepleted[i]&&totalD*i/(n-1)>d&&(connect[i-1]))
      {
         connect[i]=true;
         d=totalD*i/(n-1);
      }
   }
         std::cout<<d<<"\n";
   return d;
}
double Vtod(double V,double totalD)
{
   double impurity=-1e10*e/cm3;
   double d=sqrt(-2*epsilonGe*epsilon0*V/impurity);
   if (d>totalD)return totalD;
   //V=ax^2+c2x+c1
   //c1=0, when V=0 at x=0
   //c2=-ad, when E=dV/dx=0 at x=0, where just depleted
   //a is impurity/epsilon
   //V=-ax^2/2, solve V when x=d
   return d;

}
void UndepletedPlanarCompartoAnalytic()
{
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,500,300);
   int n=20;
   double thickess=1*cm;
   Double_t V[n], diff[n];
   for (Int_t i=1;i<n;i++) {
      V[i] = i*50;
      diff[i] = NumericalFieldGenerate(V[i],thickess)-Vtod(V[i],thickess); 
   }
   TGraph* gr = new TGraph(n,V,diff);
   gr->Draw("AC*");

}
