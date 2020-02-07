// definition of necessary units
static const double cm=1;
static const double cm3=cm*cm*cm;
static const double volt=1;
static const double C=1; // Coulomb
static const double e=1.6e-19*C; // electron charge
static const double epsilon0=8.854187817e-14*C/volt/cm; // vacuum permittivity
// https://link.springer.com/chapter/10.1007/10832182_519
static const double epsilon=15.8; // Ge dielectric constant
//______________________________________________________________________________
// (V'(r)r)'/r=a, https://www.wolframalpha.com/input/?i=(V%27(r)r)%27%2Fr%3Da
double V(double *coordinates, double *parameters)
{
   double r = coordinates[0];// there is no phi and z dependence
   double ri= parameters[0]; // inner radius
   double ro= parameters[1]; // outer radius
   double vi= parameters[2]; // inner voltage
   double vo= parameters[3]; // outer voltage
   double rho=parameters[4]; // space charge density [C/cm3]

   double a =-rho/epsilon0/epsilon;
   double c1= (vo-vi - a*(ro*ro-ri*ri)/4)/log(ro/ri);
   double c2= vo*log(ri)-vi*log(ro) - a*(ro*ro*log(ri)-ri*ri*log(ro))/4;
   c2/=log(ri)-log(ro);
   return a*r*r/4 + c1*log(r) + c2;
}
//______________________________________________________________________________
// E=-V'
double E(double *coordinates, double *parameters)
{
   double r = coordinates[0];
   double ri= parameters[0];
   double ro= parameters[1];
   double vi= parameters[2];
   double vo= parameters[3];
   double rho=parameters[4];

   double a =-rho/epsilon0/epsilon;
   double c1= (vo-vi - a*(ro*ro-ri*ri)/4)/log(ro/ri);
   return -a*r/2 - c1/r;
}
//______________________________________________________________________________
//
const int n=6; // number of curves
double rho[n]={-3.5e10*e/cm3, -1.5e10*e/cm3, 0, 
   1.5e10*e/cm3, 3.5e10*e/cm3, 6e10*e/cm3};

void drawV()
{
   TLegend *l = new TLegend(0.75,0.55,0.98,0.98);
   l->SetHeader("Impurity [cm^{-3}]");

   TF1 *fV[n]={0};
   double ri[n], ro[n], vi[n], vo[n]={0};
   for (int i=0; i<n; i++) {
      ri[i] = 0.25*cm;
      ro[i] = 1.00*cm;
      vi[i] = 2000*volt;
      fV[i] = new TF1(Form("fV%d",i), V, ri[i], ro[i], 5);
      fV[i]->SetParameters(ri[i],ro[i],vi[i],vo[i],rho[i]);
      fV[i]->SetLineStyle(i+2);
      fV[i]->SetLineColor(i+2);
      if (i+2==4) fV[i]->SetLineStyle(1);
      if (i+2==4) fV[i]->SetLineColor(1); // blue -> black
      if (i+2==5) fV[i]->SetLineColor(28); // yellow -> brown
      if (i==0) fV[i]->Draw();
      else fV[i]->Draw("same");
      // net impurity concentration = - rho/e
      l->AddEntry(fV[i],Form("%8.1e",-rho[i]/e*cm3),"l");
   }
   fV[0]->SetTitle("");
   fV[0]->GetXaxis()->SetTitle("Radial position in true coaxial detector [cm]");
   fV[0]->GetYaxis()->SetTitle("Voltage [V]");

   l->Draw();
   gPad->Print("Vrho.png");
}
//______________________________________________________________________________
//
void drawE()
{
   TCanvas *c = new TCanvas;
   TLegend *l = new TLegend(0.75,0.55,0.98,0.98);
   l->SetHeader("Impurity [cm^{-3}]");

   TF1 *fE[n]={0};
   double ri[n], ro[n], vi[n], vo[n]={0};
   for (int i=0; i<n; i++) {
      ri[i] = 0.25*cm;
      ro[i] = 1.00*cm;
      vi[i] = 2000*volt;
      fE[i] = new TF1(Form("fE%d",i), E, ri[i], ro[i], 5);
      fE[i]->SetParameters(ri[i],ro[i],vi[i],vo[i],rho[i]);
      fE[i]->SetLineStyle(i+2);
      fE[i]->SetLineColor(i+2);
      if (i+2==4) fE[i]->SetLineStyle(1);
      if (i+2==4) fE[i]->SetLineColor(1); // blue -> black
      if (i+2==5) fE[i]->SetLineColor(28); // yellow -> brown
      if (i==0) fE[i]->Draw();
      else fE[i]->Draw("same");
      // net impurity concentration = - rho/e
      l->AddEntry(fE[i],Form("%8.1e",-rho[i]/e*cm3),"l");
   }
   fE[0]->SetTitle("");
   fE[0]->GetXaxis()->SetTitle("Radial position in true coaxial detector [cm]");
   fE[0]->GetYaxis()->SetTitle("Electric field [V/cm]");

   l->Draw();
   c->Print("Erho.png");
}
//______________________________________________________________________________
//
void coaxial()
{
   gROOT->SetStyle("GeFiCa");
   gStyle->SetTitleOffset(1.2,"Y");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.01);
   gStyle->SetPadBottomMargin(0.11);

   drawV();
   drawE();
}
