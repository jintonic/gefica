// compare numerical result to analytic calculation for a 1D planar detector
using namespace GeFiCa;
void compare2analytic()
{
   Planar detector;
   detector.Bias[0]=0*volt;
   detector.Bias[1]=800*volt;
   detector.Height=1*cm;
   detector.SetAverageImpurity(1e10/cm3);

   X grid1, grid2;

   grid1.GetBoundaryConditionFrom(detector);
   grid1.SuccessiveOverRelax();

   grid2.GetBoundaryConditionFrom(detector);
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
   gn->SetTitle(";Position in a one-dimensional planar detector [cm]"
         ";Potential [V]");
   gn->GetXaxis()->SetRangeUser(0,1);
   gn->GetYaxis()->SetRangeUser(0,800);
   gn->Draw("ap");

   ga->SetLineColor(kRed);
   ga->Draw("l");

   TLegend *l = new TLegend(0.2,0.6,0.4,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR","p");
   l->Draw();

   gPad->Print("p1d.png");
}
