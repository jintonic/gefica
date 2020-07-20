using namespace GeFiCa;
// draw V, E distributions saved in an input ROOT file
void drawFields(const char *input="SPC.root")
{
   // Get data from input file
   TFile *file = new TFile(input);
   XYZ *grid = (XYZ*) file->Get("xyz");
   SquarePointContact *detector = (SquarePointContact*) file->Get("spc");
   TTree *t = grid->GetTree();
   const int n = t->GetEntries();

   t->Draw("c1:c3:v","c2>0.899&&c2<0.901","goff");
   
   TGraph2D *gv = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   gv->SetName("gv"); gv->SetNpx(100); gv->SetNpy(100); // fine bin histogram
   gv->SetTitle(";x [cm];z [cm];Voltage [V]");
   gStyle->SetNumberContours(99);
   gStyle->SetTitleOffset(1.4,"z");
   gv->Draw("surf3");

   gPad->Print("v.png");
}
