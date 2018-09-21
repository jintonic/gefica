// definition of necessary units
static const double cm=1;
static const double cm3=cm*cm*cm;
static const double volt=1;
static const double C=1; // Coulomb
static const double e=1.6e-19*C; // elementary charge
static const double epsilon0=8.854187817e-14*C/volt/cm; // vacuum permittivity
// https://link.springer.com/chapter/10.1007/10832182_519
static const double epsilon=15.8; // Ge dielectric constant
//______________________________________________________________________________
// V"(x)=a, https://www.wolframalpha.com/input/?i=V%27%27(x)%3Da
double V(double *coordinates, double *parameters)
{
   double x = coordinates[0];// there is no y and z dependence
   double x0= parameters[0]; // lower electrode
   double x1= parameters[1]; // upper electrode
   double v0= parameters[2]; // lower voltage
   double v1= parameters[3]; // upper voltage
   double rho=parameters[4]; // space charge density [C/cm3]

   double a =-rho/epsilon0/epsilon;
   double c2= (v1-v0)/(x1-x0) - a/2*(x1+x0);
   double c1= (v0*x1-v1*x0)/(x1-x0) + a/2*x0*x1;
   return a*x*x/2 + c2*x + c1;
}
//______________________________________________________________________________
// E=-V'
double E(double *coordinates, double *parameters)
{
   double x = coordinates[0];
   double x0= parameters[0];
   double x1= parameters[1];
   double v0= parameters[2];
   double v1= parameters[3];
   double rho=parameters[4];

   double a =-rho/epsilon0/epsilon;
   double c2= (v1-v0)/(x1-x0) - a/2*(x1+x0);
   return -a*x - c2;
}
//______________________________________________________________________________
//
const int n=5; // number of curves
double rho[n]={-3.5e10*e/cm3, -1.5e10*e/cm3, 0, 1.5e10*e/cm3, 3.5e10*e/cm3};

void drawV()
{
   TLegend *l = new TLegend(0.15,0.60,0.40,0.95);
   l->SetHeader("Impurity [cm^{-3}]");

   TF1 *fV[n]={0};
   double x0[n]={0}, x1[n], v0[n]={0}, v1[n];
   for (int i=0; i<n; i++) {
      x1[i] = 1*cm;
      v1[i] = 2000*volt;
      fV[i] = new TF1(Form("fV%d",i), V, x0[i], x1[i], 5);
      fV[i]->SetParameters(x0[i],x1[i],v0[i],v1[i],rho[i]);
      fV[i]->SetLineStyle(i+1);
      fV[i]->SetLineColor(i+1);
      if (i+1==5) fV[i]->SetLineColor(kMagenta); // yellow -> magenta
      if (i==0) fV[i]->Draw();
      else fV[i]->Draw("same");
      l->AddEntry(fV[i],Form("%8.1e",rho[i]/e*cm3),"l");
   }
   fV[0]->SetTitle("");
   fV[0]->GetXaxis()->SetTitle("Thickness [cm]");
   fV[0]->GetYaxis()->SetTitle("Voltage [V]");

   l->Draw();
   gPad->Print("Vx.png");
}
//______________________________________________________________________________
//
void drawE()
{
   TCanvas *c = new TCanvas;
   TLegend *l = new TLegend(0.45,0.65,0.68,0.98);
   l->SetHeader("Impurity [cm^{-3}]");

   TF1 *fE[n]={0};
   double x0[n]={0}, x1[n], v0[n]={0}, v1[n];
   for (int i=0; i<n; i++) {
      x1[i] = 1*cm;
      v1[i] = 2000*volt;
      fE[i] = new TF1(Form("fE%d",i), E, x0[i], x1[i], 5);
      fE[i]->SetParameters(x0[i],x1[i],v0[i],v1[i],rho[i]);
      fE[i]->SetLineStyle(i+1);
      fE[i]->SetLineColor(i+1);
      if (i+1==5) fE[i]->SetLineColor(kMagenta); // yellow -> magenta
      if (i==0) fE[i]->Draw();
      else fE[i]->Draw("same");
      l->AddEntry(fE[i],Form("%8.1e",rho[i]/e*cm3),"l");
   }
   fE[0]->SetTitle("");
   fE[0]->GetXaxis()->SetTitle("Thickness [cm]");
   fE[0]->GetYaxis()->SetTitle("Electric field [V/cm]");

   l->Draw();
   c->Print("Ex.png");
}
//______________________________________________________________________________
//
void planar()
{
   gROOT->SetStyle("Plain"); // pick up a good default drawing style
   // modify the default style
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetTitleOffset(1.2,"Y");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.01);
   gStyle->SetPadBottomMargin(0.11);
   
   drawV();
   drawE();
}
