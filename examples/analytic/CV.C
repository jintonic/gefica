// create C-V curve for an ideal planar detector
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
   double Cap=epsilon0*epsilonGe*1*cm*cm/d;
   return Cap*1e12;
}
//______________________________________________________________________________
//
void CV()
{
   const int n=101;
   double Cap[n],V[n];
   for (int i=0; i<n; i++) {
      V[i]=(i+1)*10*volt;
      Cap[i]=VtoC(V[i],1*cm);
   }
   // pick up a good default drawing style to modify
   gROOT->SetStyle("Plain");
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetTitleOffset(0.7,"Y");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.08);
   gStyle->SetPadTopMargin(0.01);
   // draw C VS V
   TGraph *gr  = new TGraph(n,V,Cap);
   gr->Draw("ac");
   gr->SetTitle("");
   gr->SetLineWidth(2);
   gr->GetXaxis()->SetTitle("Bias voltage [V]");
   gr->GetYaxis()->SetTitle("Capacitance [pC]");
   gPad->Print("CV.png");
}
