using namespace GeFiCa;
// search for the best relaxation factor
void optimizeRelaxationFactor()
{
   Planar detector;
   detector.Height=1*cm;
   detector.Bias[0]=0*volt;
   detector.Bias[1]=800*volt;
   detector.SetAverageImpurity(-1e10/cm3);

   gStyle->SetPadRightMargin(0.01); gStyle->SetTitleOffset(1.2,"y");
   TMultiGraph *gm = new TMultiGraph;
   TLegend *l = new TLegend(0.7,0.6,0.9,0.97);
   l->SetHeader("Grid points:");
   const int ng=6, nd=36; TGraph* g[ng]; X* grid[ng][nd];
   for (int j=0; j<ng; j++) {
      g[j] = new TGraph; g[j]->SetName(Form("g%d",ng));
      for (int i=0; i<nd; i++) {
         grid[j][i]=new X(j*100+101);
         grid[j][i]->GetBoundaryConditionFrom(detector);
         grid[j][i]->RelaxationFactor = 1.93 + i*0.002;
         grid[j][i]->MaxIterations = 10500;
         grid[j][i]->SuccessiveOverRelax();
         g[j]->SetPoint(i, grid[j][i]->RelaxationFactor,
               grid[j][i]->GetIterations());
      }
      g[j]->SetLineColor(j+1); g[j]->SetMarkerColor(j+1);
      gm->Add(g[j]);
      l->AddEntry(g[j],Form("%d",j*100+101),"pl");
   }
   gm->SetTitle(";Relaxation factor"
         ";Number of iterations for SOR to converge");
   gm->Draw("apl");
   l->Draw();
   gPad->Print("NvsRF.png");
}
