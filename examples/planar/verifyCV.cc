// Compare numerically calculated C-V curves with analytic solution
using namespace GeFiCa;
const int n=2001;
// get capacitance based on C=epsilon*A/d
double GetCfromDepletionDepth(double voltage, double thickness)
{
   // calculate fields
   Planar1D detector(n);
   detector.Thickness=thickness;
   detector.V0=voltage; // for bottom electrode
   detector.V1=0*volt; // for top electrode

   TF3 *fi = new TF3("f","-1e10"); // 1/cm3
   detector.SetImpurity(fi);
   detector.Csor=1.994;
   detector.MaxIterations=10000;
   detector.CalculatePotential(kSOR2);
   delete fi;
   
   //search for depletion depth
   double up=thickness, down=0, depth=thickness;
   while(up-down>1e-3) {
      depth=(down+up)/2;
      if(detector.GetV(depth)>0) down=depth;
      if(detector.GetV(depth)<=0) up=depth;
   }

   double A = 1*cm*1*cm;
   return epsilon*A/depth;
}
//______________________________________________________________________________
// use GefiCa::X::GetC()
double GetCfromGeFiCa(double voltage, double thickness)
{
   // calculate fields
   Planar1D detector(n);
   detector.Thickness=thickness;
   detector.V0=voltage; // for bottom electrode
   detector.V1=0*volt; // for top electrode

   TF3 *fi = new TF3("f","-1e10"); // 1/cm3
   detector.SetImpurity(fi);
   detector.Csor=1.994;
   detector.MaxIterations=10000;
   double c=detector.GetC();
   delete fi;
   
   return c;
}
//______________________________________________________________________________
//
double GetCanalytically(double voltage, double thickness)
{
   //voltage=ax^2+c2x+c1
   //c1=0, when voltage=0 at x=0
   //c2=-ad, when E=dV/dx=0 at x=0, where just depleted
   //a is rho/epsilon
   //voltage=-ax^2/2, solve voltage when x=depth
   double rho=-1e10/cm3*Qe;
   double depth=TMath::Sqrt(-2*epsilon*voltage/rho);
   if (depth>thickness) depth=thickness;

   double A = 1*cm*1*cm;
   return epsilon*A/depth;
}
//______________________________________________________________________________
//
void verifyCV()
{
   double thickness=1*cm;
   const int np=16;
   Double_t V[np], Cn[np], Ca[np],Cg[np];
   for (Int_t i=0;i<np;i++) {
      V[i] = (i+1)*50*volt; Printf("voltage: %.0f V", V[i]);
      Cn[i] = GetCfromDepletionDepth(V[i],thickness)/pF;
      Cg[i] = GetCfromGeFiCa(V[i],thickness)/pF;
      Ca[i] = GetCanalytically(V[i],thickness)/pF; 
   }

   TGraph *gn = new TGraph(np,V,Cn); gn->SetMarkerStyle(22);
   TGraph *gg = new TGraph(np,V,Cg); gg->SetMarkerStyle(24);
   TGraph *ga = new TGraph(np,V,Ca); ga->SetMarkerStyle(25);

   gStyle->SetPadLeftMargin(0.07);
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetTitleOffset(0.6,"y");
   TMultiGraph *gs = new TMultiGraph; gs->Add(gn); gs->Add(gg); gs->Add(ga);
   gs->SetTitle(";Bias [V];Capacitance [pF]"); gs->Draw("pal");

   TLegend *l = new TLegend(0.5,0.6,0.99,0.97);
   l->AddEntry(ga,"#varepsilon A/d, analytically","pl");
   l->AddEntry(gn,"#varepsilon A/d, numerically","pl");
   l->AddEntry(gg,"GeFiCa::X::GetC()","pl");
   l->Draw();
}
