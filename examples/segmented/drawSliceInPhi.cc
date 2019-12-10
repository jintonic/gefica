using namespace GeFiCa;
// Draw weighting potential of a segment of Siegfried (6 segments in phi)
void CalculateWeightingPotential()
{
   if (gSystem->Which(".","siegfried.root")) return;
   Segmented detector;
   detector.SetAverageImpurity(0);
   detector.Bias[0]=0*volt; detector.Bias[1]=1*volt; // weighting potential
   detector.Nphi=6;
   detector.SegmentId=2;

   RhoPhi grid;
   grid.SetupWith(detector);
   grid.SuccessiveOverRelax();

   TFile *output = new TFile("siegfried.root", "recreate");
   detector.Write();
   grid.Write();
   output->Close();
}
//______________________________________________________________________________
//
void DrawWeightingPotential()
{ 
   // load data
   TFile *input = new TFile("siegfried.root");
   RhoPhi *grid = (RhoPhi*) input->Get("rhophi");
   TTree *t = grid->GetTree();
   const int n = t->GetEntries();
   t->Draw("c1*cos(c2):c1*sin(c2):v","","goff");
   TGraph2D *gv = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   gv->SetName("gv"); gv->SetNpx(500); gv->SetNpy(500); // fine bin histogram

   // draw plot
   gStyle->SetPadRightMargin(0.16); gStyle->SetTitleOffset(-0.6,"z");
   TCanvas *c = new TCanvas("c", "c", 600, 520);
   c->SetLogz();
   TH2D *h = gv->GetHistogram();
   h->SetTitle(";x [cm];y [cm];weighting potential [V]");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->GetZaxis()->CenterTitle();
   h->Draw("colz");
   // widen the color palette
   c->Update(); // https://root.cern.ch/doc/master/classTPaletteAxis.html
   TPaletteAxis *palette = 
      (TPaletteAxis*) h->GetListOfFunctions()->FindObject("palette");
   palette->SetX2NDC(0.92);

   // draw segmentation scheme
   Segmented *detector = (Segmented*) input->Get("sip");
   detector->Draw("top");
   gPad->Print("sip.png");
}
//______________________________________________________________________________
//
void drawSliceInPhi()
{
   CalculateWeightingPotential();
   DrawWeightingPotential();
}
