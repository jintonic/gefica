//an example of reversed pointact detector 
{
   //basic geometry setup
   GeFiCa::ReversedCoaxialRZ *detector2 = new GeFiCa::ReversedCoaxialRZ(692,506);
   detector2->RUpperBound=3.45*GeFiCa::cm;
   detector2->RLowerBound=-3.45*GeFiCa::cm;
   detector2->ZUpperBound=5.05*GeFiCa::cm;
   detector2->PointBegin=-1.45*GeFiCa::cm;
   detector2->PointEnd=1.4500*GeFiCa::cm;
   detector2->DHole=1.0*GeFiCa::cm;
   detector2->OutterRadiusHole=1*GeFiCa::cm;
   detector2->InnerRadiusHole=0.5*GeFiCa::cm;
   detector2->removedConnorradius=1*GeFiCa::cm;
   detector2->removedConnorheight=1*GeFiCa::cm;
   TF3 *im=new TF3("f","-0.318e10+0.025e10*y");
   detector2->SetImpurity(im);
   detector2->V0=2500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;

   //calculation config
   detector2->MaxIterations=1e6;
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
   ta->Draw("c2:c1:p>>hist","","colz");
   h->GetYaxis()->SetTitle("Height [cm]");
   h->GetXaxis()->SetTitle("Radius [cm]");
  

}
