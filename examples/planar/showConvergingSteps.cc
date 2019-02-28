// Show intermediate steps of a converging process
using namespace GeFiCa;
void showConvergingSteps()
{
   Planar1D *detector = new GeFiCa::Planar1D;
   detector->SetAverageImpurity(1e10/cm3);
   detector->Thickness=1*cm;
   detector->V1=300*volt;
   detector->V0=0*volt;

   int n=3;
   TMultiGraph *mg = new TMultiGraph;
   for (int i=0;i<n;i++) {
      detector->MaxIterations=i*50;
      detector->CalculatePotential(GeFiCa::kSOR2);
      TTree *t = detector->GetTree();
      t->Draw("v:c1");
      TGraph *tg = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
      tg->SetMarkerColor(i+2);
      tg->SetMarkerStyle(6);
      mg->Add(tg);
   }
   mg->SetTitle("; Thickness [cm]; Potential [V]");

   mg->Draw("ap");
}
