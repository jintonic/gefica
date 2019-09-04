using namespace GeFiCa;
// Compare numerical result to analytic calculation
// for ideal hemispherical detector
void compare2analytic()
{
   Hemispherical detector;
   detector.Height=1*cm;
   detector.PointContactR=0.3*cm;
   detector.SetAverageImpurity(3e9/cm3);
   detector.Bias[0]=900*volt; // point contact
   detector.Bias[1]=0*volt; // outer surface

   R grid1, grid2;
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
   gn->SetTitle(";Radial position [cm];Potential [V]");
   gn->GetXaxis()->SetRangeUser(0,1);
   gn->GetYaxis()->SetRangeUser(0,900);
   gn->Draw("ap");

   ga->SetLineColor(kRed);
   ga->Draw("l");

   TLegend *l = new TLegend(0.6,0.6,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR","p");
   l->Draw();

   gPad->Print("s1d.png");
}
