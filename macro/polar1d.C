using namespace GEFICA;

void polar1d()
{
   Polar1d *detector = new Polar1d(1001);
   detector->Thickness=10*cm;
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->Create(0.3,3);
   detector->SetVoltage(2000*volt, 000*volt);
   detector->SetImpurity(1e10/cm3);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField("planar1dRK2.root");
   detector->CalculateField(EMethod::kAnalytic);
   detector->SaveField("planar1dTrue.root");
   TChain *t1=new TChain("t");
   t1->Add ("planar1dRK2.root");
   TChain *t2=new TChain("t");
   t2->Add("planar1dTrue.root");
   t1->Draw("p:c1","");
   t2->Draw("p:c1","","same");
}
