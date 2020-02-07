// Show intermediate steps of a converging process
using namespace GeFiCa;
void showConvergingSteps()
{
   Planar detector;
   detector.SetAverageImpurity(-1e10/cm3);
   detector.Height=1*cm;
   detector.Bias[0]=0*volt; // for bottom electrode
   detector.Bias[1]=300*volt; // for top electrode

   X grid; // 1D Cartesian grid
   grid.SetupWith(detector);

   // make plot
   TMultiGraph *mg = new TMultiGraph;
   TLegend *lg = new TLegend(0.66,0.135,0.99,0.57);
   const int n=5;
   int iterations[n]={0,10,20,50,164};
   for (int i=0;i<n;i++) {
      grid.MaxIterations=iterations[i];
      grid.SuccessiveOverRelax();
      TTree *t = grid.GetTree(true);
      t->Draw("v:c1","","goff");
      TGraph *tg = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
      if (i==n-1) tg->SetLineColor(7);
      else tg->SetLineColor(i+1);
      tg->SetLineStyle(i*2+1);
      tg->SetLineWidth(2);
      mg->Add(tg);
      lg->AddEntry(tg,Form("%3d iterations",iterations[i]),"l");
   }
   mg->SetTitle(";Verticle position in planar detector [cm]"
        ";Voltage [V]");
   mg->GetXaxis()->SetRangeUser(0,1);
   mg->GetYaxis()->SetRangeUser(0,350);

   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   mg->Draw("al");
   lg->Draw();

   TText *td = new TText(0.72,305,"Undepleted region");
   td->SetTextFont(132);
   td->Draw();

   TLine *tl = new TLine(0,300,1,300);
   tl->SetLineStyle(2);
   tl->Draw();

   gPad->Print("convergence.png");
}
