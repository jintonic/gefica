using namespace GEFICA;

void truee1d()
{
   TrueCoaxial1D *detector = new TrueCoaxial1D(1001);
   detector->MaxIterations=1e5;
   detector->Csor=1.9;
   detector->CalculateField(EMethod::kSOR2);
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
