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
/*
   // draw E field
   TCanvas *ce = new TCanvas;
   ce->SetLogz();
   ce->SetTopMargin(0.01);
   ce->SetRightMargin(0.11);
   ce->SetLeftMargin(0.07);

   t->Draw("c1:c2:e","","goff");
   TGraph2D *ge = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   ge->SetName("ge"); ge->SetNpx(500); ge->SetNpy(500); // fine bin histogram
   TH2D *he = ge->GetHistogram();
   he->SetTitle(";r [cm];z [cm];E [V/cm]");
   he->GetZaxis()->CenterTitle();
   he->GetZaxis()->SetTitleOffset(-0.4);
   he->GetYaxis()->SetTitleOffset(0.6);
   he->Draw("colz");
   detector->Draw();
   */
   /*
   // draw E field lines
   const int np=13;
   double x[np] = {-2.3, 0.0, 2.3, -3.4, -2.8, -0.89, 0.89, 2.8, 3.4, -3.4, 3.4, -3.4, 3.4};
   double y[np] = {5.03, 3., 5.03, 2.5, 0.03, 4.0, 4.0, 0.03, 2.5, 1.2, 1.2, 3.8,3.8};
   bool positive[np] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
   //if (detector->BoreH==0) { y[1]=5.03; y[5]=5.03; y[6]=5.03; } // no bore hole
   for (int i=0; i<np; i++)
      grid->GetFieldLineFrom(x[i]*cm, y[i]*cm, 0, positive[i])
         ->GetGraph()->Draw("p");

   ce->Print("e.png");

   // draw V field
   TCanvas *cv = new TCanvas;
   cv->SetTopMargin(0.01);
   cv->SetBottomMargin(0.05);
   cv->SetLeftMargin(0.13);
   cv->SetRightMargin(0.03);
*/
   t->Draw("c1:c3:v","c2>0.9&&c2<0.92","goff");
   
   TGraph2D *gv = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   gv->SetName("gv"); gv->SetNpx(60); gv->SetNpy(60); // fine bin histogram
   gv->SetTitle(";r [cm];z [cm];Voltage [V]");
   gStyle->SetNumberContours(99);
   gStyle->SetTitleOffset(1.4,"Z");
   gv->Draw("surf3");

  // cv->Print("v.png");
}
