using namespace GeFiCa;
// search for the best relaxation factor
void optimizeRelaxationFactor()
{
   PointContact detector; // define detector
   detector.Bias[0]=1000*kV; // bias on point contact
   detector.Bias[1]=0; // ground outer contact
   detector.SetAverageImpurity(3e9/cm3);
   detector.Radius=1*cm; detector.Height=1*cm;
   detector.PointContactR=1.0*mm; detector.PointContactH=0.1*mm;

   gStyle->SetPadRightMargin(0.01);
   TMultiGraph *gm = new TMultiGraph;
   TLegend *l = new TLegend(0.7,0.6,0.9,0.97);
   l->SetHeader("Grid points:");
   const int ng=3, nd=35; TGraph *g[ng]; RhoZ *grid[ng][nd];
   for (int j=0; j<ng; j++) {
      g[j] = new TGraph; g[j]->SetName(Form("g%d",ng));
      for (int i=0; i<nd; i++) {
         grid[j][i] = new RhoZ(j*200+50,j*200+50);
         grid[j][i]->GetBoundaryConditionFrom(detector);
         grid[j][i]->RelaxationFactor = 1.93 + i*0.002;
         grid[j][i]->MaxIterations = 8000;
         grid[j][i]->SuccessiveOverRelax();
         g[j]->SetPoint(i, grid[j][i]->RelaxationFactor,
               grid[j][i]->GetIterations());
      }
      g[j]->SetLineColor(j+1); g[j]->SetMarkerColor(j+1);
      gm->Add(g[j]);
      l->AddEntry(g[j],Form("%d#times%d",j*200+50,j*200+50),"pl");
   }
   gm->SetTitle(";Relaxation factor"
        ";Number of iterations for SOR to converge");
   gm->Draw("apl");
   l->Draw();
   gPad->Print("NvsRF4PPC.png");
}
