// definition of necessary units
static const double cm=1;
static const double cm3=cm*cm*cm;
static const double volt=1;
static const double Coulomb=1;
static const double F=Coulomb/volt; // Farat
static const double pF=1e-12*F;
static const double e=1.6e-19*Coulomb; // elementary charge
static const double epsilon0=8.854187817e-14*Coulomb/volt/cm;
// https://link.springer.com/chapter/10.1007/10832182_519
static const double epsilonGe=15.8; // Ge dielectric constant
//______________________________________________________________________________
//
double C(double V)
{
   // V"(x)=a, https://www.wolframalpha.com/input/?i=V%27%27(x)%3Da
   // => V(x) = a*x*x/2 + c2*x + c1 (a=-rho/epsilon)
   // => E(x) = -V'(x) = -a*x - c2 (a=-rho/epsilon)
   // V(x=0)=0 => c1=0
   // E(x=d)=0 => c2=-ad (d is the thickness of the depleted region)
   // Hence V(d) = -ad^2/2, hence d = sqrt(-2V/a)
   double rho=1e10*e/cm3; // space charge density
   double a=-rho/epsilon0/epsilonGe;
   double d=sqrt(-2*V/a);
   double thickness=1*cm; // detector thickness
   if (d>thickness) d=thickness;
   double area=1*cm*1*cm; // unit area
   return epsilon0*epsilonGe*area/d; // C = epsilon * area / d
}
//______________________________________________________________________________
//
void CV()
{
   const int n=101; double c[n],v[n];
   for (int i=0; i<n; i++) {
      v[i]=(i+1)*10;
      c[i]=C(v[i]*volt)/pF;
   }
   TGraph *gr = new TGraph(n,v,c);

   gROOT->SetStyle("GeFiCa");
   gStyle->SetTitleOffset(0.7,"Y");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.08);
   gStyle->SetPadTopMargin(0.01);

   // draw C VS V
   gr->Draw("ac");
   gr->SetTitle(";Bias voltage [V];Capacitance [pC]");
   gr->SetLineWidth(2);

   gPad->Print("CV.png");
}
