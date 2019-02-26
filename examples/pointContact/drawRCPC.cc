/**
 * \file drawRCPC.cc
 * \example pointContact/drawRCPC.cc
 * \brief draw fields of a 2D reversed coaxial point contact detector
 */
{
   //basic geometry setup
   GeFiCa::ReversedCoaxialRZ *detector2 = new GeFiCa::ReversedCoaxialRZ(692,506);
   detector2->Radius=3.45*GeFiCa::cm;
   detector2->Z=5.05*GeFiCa::cm;
   detector2->PointContactR=1.45*GeFiCa::cm;
   detector2->HoleZ=2.0*GeFiCa::cm;
   detector2->HoleOutterR=1.2*GeFiCa::cm;
   detector2->HoleInnerR=0.5*GeFiCa::cm;
   detector2->ConnorZ=1.2*GeFiCa::cm;
   detector2->ConnorLength=1.3*GeFiCa::cm;
   TF3 *im=new TF3("f","-0.318e10+0.025e10*y");
   detector2->SetImpurity(im);
   detector2->V0=2500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;

   //calculation config
   detector2->MaxIterations=2e3;
   detector2->Precision=1e-8;
   detector2->Csor=1.992;

   //calculate field
   detector2->CalculatePotential(GeFiCa::kSOR2);
   detector2->SaveField("rcpcSOR2.root");
   
   //prepare drawing style
   gStyle->SetOptTitle(kTRUE);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadRightMargin(1.01);
   gStyle->SetPadLeftMargin(0.0999999999);
   gStyle->SetLabelFont(22,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetTitleFont(22,"XY");
   gStyle->SetLegendFont(22);

   //generate graphics
   TH2F* h= new TH2F("hist","",10,-3.45,3.45,10,0,5.05);
   TChain *ta = new TChain("t");
   ta->Add("rcpcSOR2.root");
   ta->Draw("c2:c1:v>>hist","","colz");
   h->GetYaxis()->SetTitle("Height [cm]");
   h->GetXaxis()->SetTitle("Radius [cm]");
}