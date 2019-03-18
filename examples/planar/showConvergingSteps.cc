// Show intermediate steps of a converging process
using namespace GeFiCa;
void showConvergingSteps()
{
   // setup grid
   Planar1D *detector = new GeFiCa::Planar1D;
   detector->SetAverageImpurity(1e10/cm3);
   detector->Thickness=1*cm;
   detector->V1=300*volt;
   detector->V0=0*volt;

   // prepare drawing style
   gROOT->SetStyle("Plain"); // pick up a good default style to modify
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.02);

   // make plot
   TMultiGraph *mg = new TMultiGraph;
   TLegend *lg = new TLegend(0.66,0.135,0.99,0.57);
   const int n=5;
   int nIter[n]={0,10,20,50,234};
   for (int i=0;i<n;i++) {
      detector->MaxIterations=nIter[i];
      detector->SuccessiveOverRelax();
      TTree *t = detector->GetTree(true);
      t->Draw("v:c1","","goff");
      TGraph *tg = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
      if (i==n-1) tg->SetLineColor(7);
      else tg->SetLineColor(i+1);
      tg->SetLineStyle(i*2+1);
      tg->SetLineWidth(2);
      mg->Add(tg);
      lg->AddEntry(tg,Form("%3d iterations",nIter[i]),"l");
   }
   mg->SetTitle(";Position in a one dimensional planar detector [cm]"
        ";Potential [V]");
   mg->GetXaxis()->SetRangeUser(0,1);
   mg->GetYaxis()->SetRangeUser(0,350);
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
