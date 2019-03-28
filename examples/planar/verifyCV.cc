// Compare numerically calculated C-V curves with analytic solution
using namespace GeFiCa;
const int n=401;
// get capacitance based on C=epsilon*A/d
double GetCfromDepletionDepth(double voltage, double height)
{
   Planar detector;
   detector.Height=height;
   detector.Bias[0]=voltage; // for bottom electrode
   detector.Bias[1]=0*volt; // for top electrode
   detector.SetAverageImpurity(-1e10/cm3);

   X grid(n);
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.986;
   grid.SuccessiveOverRelax();
   
   //search for depletion depth
   double up=height, down=0, depth=height;
   while(up-down>1e-3) {
      depth=(down+up)/2;
      gDebug=1;
      if(grid.GetV(depth)==0) down=depth; else up=depth;
      gDebug=0;
   }

   double A = 1*cm*1*cm;
   return epsilon*A/depth;
}
//______________________________________________________________________________
// use GefiCa::X::GetC()
double GetCfromGeFiCa(double voltage, double height)
{
   // calculate fields
   Planar detector;
   detector.Height=height;
   detector.Bias[0]=voltage; // for bottom electrode
   detector.Bias[1]=0*volt; // for top electrode
   detector.SetAverageImpurity(-1e10/cm3);

   X grid(n);
   grid.GetBoundaryConditionFrom(detector);
   grid.RelaxationFactor=1.986;
   return grid.GetC();
}
//______________________________________________________________________________
//
double GetCanalytically(double voltage, double height)
{
   //voltage=ax^2+c2x+c1
   //c1=0, when voltage=0 at x=0
   //c2=-ad, when E=dV/dx=0 at x=0, where just depleted
   //a is rho/epsilon
   //voltage=-ax^2/2, solve voltage when x=depth
   double rho=1e10/cm3*Qe;
   double depth=TMath::Sqrt(2*epsilon*voltage/rho);
   if (depth>height) depth=height;

   double A = 1*cm*1*cm;
   return epsilon*A/depth;
}
//______________________________________________________________________________
//
void verifyCV()
{
   double height=1*cm;
   const int np=16;
   Double_t V[np], Cn[np], Ca[np],Cg[np];
   for (Int_t i=0;i<np;i++) {
      V[i] = (i+2)*50*volt; Printf("voltage: %.0f V", V[i]);
      Cn[i] = GetCfromDepletionDepth(V[i],height)/pF;
      Cg[i] = GetCfromGeFiCa(V[i],height)/pF;
      Ca[i] = GetCanalytically(V[i],height)/pF; 
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
