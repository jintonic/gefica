#include <TFile.h>
#include <Planar1D.h>
using namespace GEFICA;

void planar(const char* fin="input.root")
{
   detector->field = new Field2D(10);
   field->Load(fin);

   Planar *detector = new Planar;
   detector->voltage = -2000;
   detector->impurity = new TF1("impurity","3*x+1",0,10);
   detector->height = 1;

   detector->SetNumericalMethod("asjdfkas");
   detector->UpdateField();
   field->Write("output.root");

   TFile *out = new TFile("output.root","recreate");
   TTree *t = new TTree("t","weighting potential");
   t->Branch("x",x,"x[n1]/D");

   for (unsigned short i=0; i<n1; i++) {
      t->Fill();
   }
}
