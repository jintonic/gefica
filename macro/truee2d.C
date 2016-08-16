using namespace GEFICA;

void truee2d()
{
   TrueCoaxial2D *detector = new TrueCoaxial2D(101,101);
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->CalculateField(EMethod::kSOR2);
   detector->SaveField("planar1dRK2.root");

   TChain *t1=new TChain("t");
   t1->Add ("planar1dRK2.root");
   TChain *t2=new TChain("t");
   t1->Draw("p:c1:c2");
}
