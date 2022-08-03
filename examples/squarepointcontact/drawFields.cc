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

   gStyle->SetNumberContours(99);

   // draw E field
   TCanvas *ce = new TCanvas;
   //ce->SetLogz();
   ce->SetTopMargin(0.02);
   ce->SetRightMargin(0.12);
   ce->SetLeftMargin(0.08);

   t->Draw("c1:c3:e","c2>0.899&&c2<0.901","goff");
   TGraph2D *ge = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   ge->SetName("ge"); ge->SetNpx(200); ge->SetNpy(200); // fine bin histogram
   TH2D *he = ge->GetHistogram();
   he->SetTitle(";x [cm];z [cm];E [V/cm]");
   he->GetZaxis()->CenterTitle();
   he->GetZaxis()->SetTitleOffset(-0.4);
   he->GetYaxis()->SetTitleOffset(0.75);
   he->Draw("colz");
   ce->Print("e.png");

   // draw V field
   TCanvas *cv = new TCanvas;
   cv->SetTopMargin(0.01);
   cv->SetBottomMargin(0.05);
   cv->SetLeftMargin(0.13);
   cv->SetRightMargin(0.03);
   
   t->Draw("c1:c3:v","c2>0.899&&c2<0.901","goff");
   TGraph2D *gv = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   gv->SetName("gv"); gv->SetNpx(100); gv->SetNpy(100); // fine bin histogram
   gv->SetTitle(";x [cm];z [cm];Voltage [V]");
   gv->Draw("surf3");

   gPad->Print("v.png");
}
