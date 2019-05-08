// Compare numerically calculated C-V curves with analytic solutions
using namespace GeFiCa;
const int n=401;
// get capacitance based on C=2pi*epsilon/ln(radius/boreR), ref. to
// http://hyperphysics.phy-astr.gsu.edu/hbase/electric/capcyl.html
double GetCfromDepletionRadius(double voltage, double radius, double boreR)
{
   TrueCoaxial detector;
   detector.Radius=radius;
   detector.BoreR=boreR;
   detector.Bias[0]=0*volt; // for inner electrode
   detector.Bias[1]=voltage; // for outer electrode
   detector.SetAverageImpurity(-2e9/cm3); // deplete from inside out

   Rho grid(n);
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.986;
   grid.Precision=1e-8;
   grid.SuccessiveOverRelax();
   
   //search for depletion radius
   double up=radius, down=boreR, depletionR=radius;
   while(up-down>1e-7*mm) {
      depletionR=(down+up)/2;
      if (grid.GetV(depletionR)<voltage) down=depletionR; else up=depletionR;
   }

   return epsilon*3.14159*2/log(depletionR/boreR); // c per unit length in cm
}
//______________________________________________________________________________
// use GefiCa::X::GetC()
double GetCfromGeFiCa(double voltage, double radius,double boreR)
{
   // calculate fields
   TrueCoaxial detector;
   detector.Radius=radius;
   detector.BoreR=boreR;
   detector.Bias[0]=0*volt; // for bottom electrode
   detector.Bias[1]=voltage; // for top electrode
   detector.SetAverageImpurity(-2e9/cm3);

   Rho grid(n);
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.986;
   grid.Precision=1e-8;

   return grid.GetC();
}
//______________________________________________________________________________
// search for depleted radius recursively
double FindR(double voltage,double upR,double downR,double boreR,double rho)
{
   if (upR<downR+1e-7*mm) return upR; // find the required R
   double r=(upR+downR)/2;
   //voltage= -rho*r*r/4 + c1*log(r) + c2, where rho=-2e9/cm3*Qe/epsilon
   //c1=-rho*r^2/2 <= when E(r)=0, the detector is depleted from boreR to r
   //c2=rho*boreR*boreR/4 - c1*log(boreR) <= V(boreR)=0
   double c1=-rho*r*r/2;
   double c2=rho*boreR*boreR/4-c1*log(boreR);
   double V=rho*r*r/4+c1*log(r)+c2;
   if (V>voltage) return FindR(voltage,r,downR,boreR,rho);
   else return FindR(voltage,upR,r,boreR,rho);
}
//______________________________________________________________________________
//
double GetCanalytically(double voltage, double radius,double boreR)
{
   double rho=-2e9/cm3*Qe/epsilon;
   double depth=FindR(voltage,radius,boreR,boreR,rho);
   if (depth>radius) depth=radius;
   return epsilon*3.14159*2/log(depth/boreR); // c per unit length in cm
}
//______________________________________________________________________________
//
void verifyCV()
{
   double radius=3.5*cm;
   double boreR=0.5*cm;
   const int np=25;
   Double_t V[np], Cn[np], Ca[np],Cg[np];
   for (Int_t i=0;i<np;i++) {
      V[i] = (i+2)*100*volt; Printf("voltage: %.0f V", V[i]);
      Cn[i] = GetCfromDepletionRadius(V[i],radius,boreR)/pF;
      Cg[i] = GetCfromGeFiCa(V[i],radius,boreR)/pF;
      Ca[i] = GetCanalytically(V[i],radius,boreR)/pF; 
   }

   TGraph *gn = new TGraph(np,V,Cn); gn->SetMarkerStyle(22);
   TGraph *gg = new TGraph(np,V,Cg); gg->SetMarkerStyle(24);
   TGraph *ga = new TGraph(np,V,Ca); ga->SetMarkerStyle(25);

   gStyle->SetPadLeftMargin(0.09);
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetTitleOffset(0.85,"y");
   TMultiGraph *gs = new TMultiGraph; gs->Add(gn); gs->Add(gg); gs->Add(ga);
   gs->SetTitle(";Bias [V];Capacitance per unit length [pF/cm]");
   gs->Draw("pac");

   TLegend *l = new TLegend(0.4,0.7,0.97,0.97);
   l->AddEntry(ga,"2#pi#varepsilon/ln(R_{out}/R_{in}), analytically","pl");
   l->AddEntry(gn,"2#pi#varepsilon/ln(R_{out}/R_{in}), numerically","pl");
   l->AddEntry(gg,"GeFiCa::X::GetC()","pl");
   l->Draw();
}
