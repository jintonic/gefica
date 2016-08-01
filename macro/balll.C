using namespace GEFICA;

void balll()
{
   ball *detector = new ball(11,11,11);
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->Create(0.3,3);
   detector->SetVoltage(0*volt,10*volt);
   detector->CalculateField(EMethod::kRK2);
   detector->SaveField("planar1dRK2.root");

   Polar1d *detector2 = new Polar1d(11);
   detector2->Thickness=10*cm;
   detector2->MaxIterations=1e5;
   detector2->Csor=1.9;
   detector2->Create(0.3,3);
   detector2->SetVoltage(0*volt, 10*volt);
   detector2->CalculateField(EMethod::kAnalytic);
   detector2->SaveField("planar1dTrue.root");

   TChain *t1=new TChain("t");
   t1->Add ("planar1dRK2.root");
   TChain *t2=new TChain("t");
   t2->Add("planar1dTrue.root");
   t1->Draw("p:c1","c1<1e3");
   t2->Draw("p:c1","","same");
}
