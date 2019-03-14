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
   gStyle->SetTitleOffset(0.75,"Y");
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadLeftMargin(0.08);

   // generate plots
   t->Draw("c1:c2:log(v)","","goff");
   TGraph2D *gv = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   gv->SetName("gv"); gv->SetNpx(500); gv->SetNpy(500); // fine bin histogram
   TH2D *hv = gv->GetHistogram();
   hv->SetTitle(";Radius [cm];Height [cm];Potential [V]");
   hv->GetZaxis()->CenterTitle();
   hv->Draw("colz");

   // draw E field lines
   const int np=12;
   TGraph *gl[np];
   double x[np] = {-25, 0.01, 25, -4, -3, -2, -1, 0.01, 1, 2, 3, 4};
   double y[np] = {49.9, 24.9, 49.9, 2, 2, 2, 2, 2, 2, 2, 2, 2};
   bool positive[np] = {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
   for (int i=0; i<np; i++) {
      gDebug=1;
      gl[i] = detector->GetFieldLineFrom(x[i]*mm, y[i]*mm, positive[i]);
      gDebug=0;
      gl[i]->Draw("l");
   }

   TCanvas *ce = new TCanvas;
   ce->SetLogz();
   t->Draw("c1:c2:e","","goff");
   TGraph2D *ge = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   ge->SetName("ge"); ge->SetNpx(500); ge->SetNpy(500); // fine bin histogram
   TH2D *he = ge->GetHistogram();
   he->SetTitle(";Radius [cm];Height [cm];E [V/cm]");
   he->GetZaxis()->CenterTitle();
   he->Draw("colz");
   for (int i=0; i<np; i++) gl[i]->Draw("p");
}
