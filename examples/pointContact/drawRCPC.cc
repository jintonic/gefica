// demonstrate the advantage of 2D reversed coaxial point contact detector
using namespace GeFiCa;
void drawRCPC()
{
   const int nr=500, nz=251;
   PointContactDZ *rcpc = new PointContactDZ(nr,nz);
   rcpc->Height=1.00*cm;
   rcpc->Radius=1.00*cm;
   rcpc->PointContactH=0.1*mm;
   rcpc->PointContactR=1.0*mm;
   rcpc->HoleH=0.5*cm;
   rcpc->HoleInnerR=0.2*cm;
   rcpc->HoleOuterR=0.3*cm;
   rcpc->TaperH=2*mm;
   rcpc->TaperW=2*mm;
   rcpc->CornerH=2*mm;
   rcpc->CornerW=2*mm;
   rcpc->WrapArroundR=0.5*cm;
   TF3 *im = new TF3("f","0.318e10-0.025e10*y");
   rcpc->SetImpurity(im);
   rcpc->V0=0*volt;
   rcpc->V1=1000*volt;

   //calculate field
   rcpc->Precision=1e-8;
   rcpc->Csor=1.992;
   rcpc->CalculatePotential(kSOR2);
   
   // prepare drawing style
   gROOT->SetStyle("Plain"); // pick up a good drawing style to modify
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XYZ");
   gStyle->SetTitleFont(132,"XYZ");
   gStyle->SetLabelSize(0.05,"XYZ");
   gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(0.8,"Y");
   gStyle->SetTitleOffset(-0.4,"Z");
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadLeftMargin(0.08);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetOptStat(0);

   //generate plot
   TTree *t = rcpc->GetTree();
   t->Draw("c1:c2:sqrt(e1*e1+e2*e2)","","goff"); //it creates TH3 instead of TH2
   const int n=t->GetSelectedRows();
   TGraph2D *g = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
   g->SetNpx(nr-1); g->SetNpy(nz-1); // set numbers of bins for TH2
   TH2D *h = g->GetHistogram();
   h->SetTitle(";Radius [cm];Height [cm];E [V/cm]");
   h->GetZaxis()->CenterTitle(); // only works for h not g
   h->Draw("colz"); // "colz" only works for TH2
   gPad->SetLogz();
}
