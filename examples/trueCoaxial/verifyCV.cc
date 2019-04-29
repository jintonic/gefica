// Compare numerically calculated C-V curves with analytic solution

using namespace GeFiCa;
const int n=401;
// get capacitance based on C=epsilon*A/d
double GetCfromDepletionDepth(double voltage, double radius,double borer)
{
   TrueCoaxial detector;
   detector.Radius=radius;
   detector.BoreR=borer;
   detector.Bias[0]=0*volt; // for bottom electrode
   detector.Bias[1]=voltage; // for top electrode
   detector.SetAverageImpurity(-1e10/cm3); // deplete from bottom

   Rho grid(n);
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.986;
   grid.SuccessiveOverRelax();
   
   //search for depletion depth
   double up=radius, down=borer, depth=radius;
   while(up-down>1e-3) {
      depth=(down+up)/2;
      if (grid.GetV(depth)<voltage) down=depth; else up=depth;
   }

   double A = 1*cm*1*cm;
   return epsilon*A/depth;
}
//______________________________________________________________________________
// use GefiCa::X::GetC()
double GetCfromGeFiCa(double voltage, double radius,double borer)
{
   // calculate fields
   TrueCoaxial detector;
   detector.Radius=radius;
   detector.BoreR=borer;
   detector.Bias[0]=0*volt; // for bottom electrode
   detector.Bias[1]=voltage; // for top electrode
   detector.SetAverageImpurity(-1e10/cm3);

   Rho grid(n);
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.986;
   return grid.GetC();
}
//______________________________________________________________________________
//
double findr(double voltage,double topr,double downr, double rho)
{
   if(topr<downr)return topr;
   double r=(topr+downr)/2;
   double c1=-rho*r*r/2;
   double V=1/4*rho*r*r+c1*log(r);
   cout<<V<<endl;
   if(V>voltage)return findr(voltage,r-1e-5,downr,rho);
   else if(V<voltage)return findr(voltage,topr,r+1e-5,rho);
   else return topr;


}
//______________________________________________________________________________
//
double GetCanalytically(double voltage, double radius,double borer)
{
   //voltage=ax^2+c2x+c1
   //c1=0, when voltage=0 at x=0
   //c2=-ad, when E=dV/dx=0 at x=0, where just depleted
   //a is rho/epsilon
   //voltage=-ax^2/2, solve voltage when x=depth
   double rho=-1e10/cm3*Qe;
   double depth=findr(voltage,radius,borer,rho);
   if (depth>radius) depth=radius;

   double L = 1*cm;
   return epsilon*L*3.14159*2/log(depth/borer);
}
//______________________________________________________________________________
//
void verifyCV()
{
   double radius=1*cm;
   double borer=0.1*cm;
   const int np=16;
   Double_t V[np], Cn[np], Ca[np],Cg[np];
   for (Int_t i=0;i<np;i++) {
      V[i] = (i+2)*50*volt; Printf("voltage: %.0f V", V[i]);
      Cn[i] = GetCfromDepletionDepth(V[i],radius,borer)/pF;
      Cg[i] = GetCfromGeFiCa(V[i],radius,borer)/pF;
      Ca[i] = GetCanalytically(V[i],radius,borer)/pF; 
   }

   TGraph *gn = new TGraph(np,V,Cn); gn->SetMarkerStyle(22);
   TGraph *gg = new TGraph(np,V,Cg); gg->SetMarkerStyle(24);
   TGraph *ga = new TGraph(np,V,Ca); ga->SetMarkerStyle(25);

   gStyle->SetPadLeftMargin(0.09);
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetTitleOffset(0.85,"y");
   TMultiGraph *gs = new TMultiGraph; gs->Add(gn); gs->Add(gg); gs->Add(ga);
   gs->SetTitle(";Bias [V];Capacitance per cm^{2} [pF]"); gs->Draw("pac");

   TLegend *l = new TLegend(0.5,0.6,0.97,0.97);
   l->AddEntry(ga,"#varepsilon A/d, analytically","pl");
   l->AddEntry(gn,"#varepsilon A/d, numerically","pl");
   l->AddEntry(gg,"GeFiCa::X::GetC()","pl");
   l->Draw();
}
