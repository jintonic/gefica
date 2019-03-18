using namespace GeFiCa;
// Draw weighting potential of a segment of Siegfried (6 segments in phi)
void CalculateWeightingPotential()
{
   if (gSystem->Which(".","siegfried.root")) return;
   SegmentedInPhi *siegfried = new SegmentedInPhi;
   siegfried->V0=0*volt;
   siegfried->V1=1*volt;
   siegfried->Nseg=6;
   siegfried->SegmentId=2;
   siegfried->SuccessiveOverRelax();
   TFile *output = new TFile("siegfried.root", "recreate");
   siegfried->Write();
   output->Close();
}
//______________________________________________________________________________
//
void DrawWeightingPotential()
{ 
   // load data
   TFile *input = new TFile("siegfried.root","update");
   SegmentedInPhi *siegfried = (SegmentedInPhi*) input->Get("sip");
   TTree *t = siegfried->GetTree();
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
   double r=siegfried->OuterR, ri=siegfried->InnerR;
   double x=r*cos(3.14/3), y=r*sin(3.14/3);

   TLine *l1 = new TLine(-r,0,r,0);
   l1->SetLineColor(kBlack); l1->SetLineStyle(kDashed); l1->Draw();
   TLine *l2 = new TLine(-x,-y,x,y);
   l2->SetLineColor(kBlack); l2->SetLineStyle(kDashed); l2->Draw();
   TLine *l3 = new TLine(-x,y,x,-y);
   l3->SetLineColor(kBlack); l3->SetLineStyle(kDashed); l3->Draw();

   TEllipse *e1 = new TEllipse(0,0,ri,ri); e1->Draw();
   TEllipse *e2 = new TEllipse(0,0,r,r); e2->SetFillStyle(0); e2->Draw();

   gPad->Print("sip.png");
}
//______________________________________________________________________________
//
void drawSliceInPhi()
{
   CalculateWeightingPotential();
   DrawWeightingPotential();
}
