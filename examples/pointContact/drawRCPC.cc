// demonstrate the advantage of 2D reversed coaxial point contact detector
using namespace GeFiCa;
void drawRCPC()
{
   //const int nr=692, nz=506;
   const int nr=300, nz=301;
   PointContactDZ *rcpc = new PointContactDZ(nr,nz);
   rcpc->Height=5.05*cm;
   rcpc->Radius=3.45*cm;
   rcpc->PointContactH=0.1*mm;
   rcpc->PointContactR=1.0*mm;
   rcpc->HoleH=2.0*cm;
   rcpc->HoleInnerR=0.5*cm;
   rcpc->HoleOuterR=1.2*cm;
   rcpc->TaperH=4*mm;
   rcpc->TaperW=4*mm;
   rcpc->CornerH=4*mm;
   rcpc->CornerW=4*mm;
   rcpc->V1=0*volt;
   rcpc->V0=3000*volt;
   TF3 *im = new TF3("f","-0.318e10+0.025e10*y");
   rcpc->SetImpurity(im);

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
   gStyle->SetTitleOffset(0.6,"Y");
   gStyle->SetTitleOffset(-0.4,"Z");
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadLeftMargin(0.07);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetOptStat(0);

   //generate graphics
   TTree *t = rcpc->GetTree();
   TH2D *h = new TH2D("h", ";Radius [cm];Height [cm];Potential [V]",
         nr-1, -rcpc->Radius, rcpc->Radius, nz-1, 0, rcpc->Height);
   h->GetZaxis()->CenterTitle();
   t->Draw("c2:c1:v>>h","","colz");
}
