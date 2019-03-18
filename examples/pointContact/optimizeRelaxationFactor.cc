using namespace GeFiCa;
// search for the best Csor
void optimizeCsor()
{
   gStyle->SetPadRightMargin(0.01);
   TMultiGraph *gm = new TMultiGraph;
   TLegend *l = new TLegend(0.65,0.6,0.85,0.97);
   l->SetHeader("Grid points:");
   const int ng=6, nd=34; TGraph *g[ng]; PointContactDZ *detector[ng][nd];
   for (int j=0; j<ng; j++) {
      g[j] = new TGraph; g[j]->SetName(Form("g%d",ng));
      for (int i=0; i<nd; i++) {
         detector[j][i] = new PointContactDZ(j*100+50,j*100+50);
         detector[j][i]->V0=1000*volt;
         detector[j][i]->V1=0*volt;
         detector[j][i]->SetAverageImpurity(3e9/cm3);

         detector[j][i]->Csor = 1.93 + i*0.002;
         detector[j][i]->MaxIterations = 6000;
         detector[j][i]->CalculatePotential();
         detector[j][i]->Gsor->SetName(Form("gsor%d",i));
         g[j]->SetPoint(i, detector[j][i]->Csor, detector[j][i]->GetNsor());
      }
      g[j]->SetLineColor(j+1); g[j]->SetMarkerColor(j+1);
      gm->Add(g[j]);
      l->AddEntry(g[j],Form("%d#times%d",j*100+50,j*100+50),"pl");
   }
   gm->SetTitle(";Constant of SOR (Csor)"
        ";Number of iterations for SOR to converge");
   gm->Draw("apl");
   l->Draw();
   gPad->Print("NvsCsor.png");
}
