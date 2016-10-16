using namespace GeFiCa;

void planar2d()
{
   Planar2D *detector = new Planar2D(101,101);
   detector->MaxIterations=1e5;
   detector->Csor=1;
   detector->CalculateField(EMethod::kSOR2);
   detector->SaveField("planar1dRK2.root");

   TCanvas *c = new TCanvas;
   TChain *t1 = new TChain("t");
   t1->Add("planar1dRK2.root");
   t1->Draw("p:c1:c2");
   return c;


}
