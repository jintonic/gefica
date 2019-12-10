using namespace GeFiCa;
// compare numerical result to analytic calculation for a true coaxial detector
void compare2analytic()
{
   TrueCoaxial detector;
   detector.Bias[0]=3000*volt; // core
   detector.Bias[1]=0*volt; // outer surface electrode
   detector.BoreR=0.5*cm;
   detector.Radius=3.0*cm;
   detector.SetAverageImpurity(-3e9/cm3); // n-type

   Rho grid1, grid2;
   grid1.SetupWith(detector);
   grid1.SuccessiveOverRelax();
   grid2.SetupWith(detector);
   grid2.SolveAnalytically();

   // generate graphics
   TTree *tn = grid1.GetTree();
   tn->Draw("v:c1","","goff");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());
   TTree *ta = grid2.GetTree();
   ta->Draw("v:c1","","goff");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // compare numerical result to analytic calculation
   gStyle->SetPadLeftMargin(0.11); gStyle->SetPadRightMargin(0.01);
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Radial position [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

   TLegend *l = new TLegend(0.5,0.65,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR","p");
   l->Draw();

   gPad->Print("tc.png");
}
