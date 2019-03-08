using namespace GeFiCa;
// draw V, E distributions saved in an input ROOT file
void drawFields(const char *input="ppc.root")
{
   // Get data from input file
   TFile *file = new TFile(input, "update"); // "update" allows creation of TH3
   PointContactDZ *detector = (PointContactDZ*) file->Get("pcdz");
   TTree *t = detector->GetTree(true);
   const int n = t->GetEntries();

   // fine tune margins, etc.
   gStyle->SetTitleOffset(0.6,"Y");
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadLeftMargin(0.07);

   // generate plots
   TCanvas *cv = new TCanvas;
   //cv->SetLogz();
   t->Draw("c1:c2:v","","goff");
   TGraph2D *gv = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   gv->SetName("gv"); gv->SetNpx(500); gv->SetNpy(500); // fine bin histogram
   TH2D *hv = gv->GetHistogram();
   hv->SetTitle(";Radius [cm];Height [cm];Potential [V]");
   hv->GetZaxis()->CenterTitle();
   hv->Draw("colz");

   TCanvas *ce = new TCanvas;
   ce->SetLogz();
   t->Draw("c1:c2:e","","goff");
   TGraph2D *ge = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   ge->SetName("ge"); ge->SetNpx(500); ge->SetNpy(500); // fine bin histogram
   TH2D *he = ge->GetHistogram();
   he->SetTitle(";Radius [cm];Height [cm];E [V/cm]");
   he->GetZaxis()->CenterTitle();
   he->Draw("colz");

   gDebug=1;
   TGraph *g1 = detector->GetFieldLineFrom(3*cm, 1*cm);
   g1->Draw("pl");
}
