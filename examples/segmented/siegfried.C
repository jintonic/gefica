// Draw weighting potential of a segment of Siegfried (6 segments in phi)
void CalculateWeightingPotential()
{
   GeFiCa::Segmented2D *siegfried = new GeFiCa::Segmented2D;
   siegfried->V0=0*GeFiCa::volt;
   siegfried->V1=1*GeFiCa::volt;
   siegfried->SegmentNum=6;
   siegfried->CalculatePotential(GeFiCa::kSOR2,3);
   siegfried->SaveField("siegfried.root");
}
//______________________________________________________________________________
//
void DrawWeightingPotential()
{ 
   // pick up a good default drawing style to modify
   gROOT->SetStyle("Plain");
   gStyle->SetTitleFont(132,"XYZ"); gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetLabelFont(132,"XYZ"); gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadLeftMargin(0.08);
   gStyle->SetPadTopMargin(0.02);
   // create a smoother palette than the default one
   const int nRGBs = 5;
   const int nCont = 255;
   double stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   double red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   double green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   double blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
   gStyle->SetNumberContours(nCont);
   // make a square canvas
   TCanvas *c = new TCanvas("c", "c", 450, 450);
   c->SetLogz();
   // load data and draw
   TChain *t = new TChain("t");
   t->Add("siegfried.root");
   t->Draw("c1*sin(c2):c1*cos(c2):p","","colz");
   // fine tune plot
   TH2F *h = (TH2F*) gPad->GetPrimitive("htemp");
   h->SetTitle(";x [cm];y [cm];weighting potential [V]");
   h->GetYaxis()->SetTitleOffset(0.7);
   h->GetZaxis()->SetTitleOffset(-0.4);
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->GetZaxis()->CenterTitle();
   // widen the color palette
   gPad->Update(); // create the palette by forcely drawing the plot
   TPaletteAxis *palette = 
      (TPaletteAxis*) h->GetListOfFunctions()->FindObject("palette");
   palette->SetX2NDC(0.94);
   h->Draw("colz"); // let the new setup take effect
   // draw segmentation scheme
   GeFiCa::Segmented2D *siegfried = new GeFiCa::Segmented2D();
   double r=siegfried->RUpperBound, x=r*cos(3.14/3), y=r*sin(3.14/3);
   TLine *l1 = new TLine(-r,0,r,0);
   l1->SetLineColor(kWhite);
   l1->SetLineStyle(kDashed);
   l1->Draw();
   TLine *l2 = new TLine(-x,-y,x,y);
   l2->SetLineColor(kWhite);
   l2->SetLineStyle(kDashed);
   l2->Draw();
   TLine *l3 = new TLine(-x,y,x,-y);
   l3->SetLineColor(kWhite);
   l3->SetLineStyle(kDashed);
   l3->Draw();
}
//______________________________________________________________________________
//
void siegfried()
{
   if (1||gSystem->Which(".","siegfried.root")==0) CalculateWeightingPotential();
   DrawWeightingPotential();
}
