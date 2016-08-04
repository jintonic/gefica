using namespace GEFICA;

TCanvas* planar1d()
{
   Planar1D *detector = new Planar1D(101);
   detector->Thickness=10*cm;
   detector->MaxIterations=1e5;
   detector->Csor=1;
   detector->SetVoltage(0*volt, 2000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField("planar1dRK2.root");
   detector->CalculateField(EMethod::kAnalytic);
   detector->SaveField("planar1dTrue.root");

   TCanvas *c = new TCanvas;
   TChain *t1 = new TChain("t");
   t1->Add("planar1dRK2.root");
   t1->Draw("p:c1");
   return c;
}
