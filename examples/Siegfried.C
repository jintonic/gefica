
{
   GeFiCa::Siegfried2D *wp = new GeFiCa::Siegfried2D(301,901);
   wp->MaxIterations=1e5;
   wp->Csor=1.999;
   wp->V0=0*GeFiCa::volt;
   wp->V1=1*GeFiCa::volt;
   //wp->Impurity="1e10";
   //wp->CalculatePotential(GeFiCa::kSOR2);
   //wp->SaveField("Siegfried2D.root");
   const Int_t NRGBs = 5;
   const Int_t NCont = 255;

   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
   TCanvas * c1 = new TCanvas("c", "c", 450, 450);
   TChain *tn = new TChain("t");
   tn->Add("Siegfried2D.root");
   tn->Draw("c1*sin(c2):c1*cos(c2):p","","colz ");

}
