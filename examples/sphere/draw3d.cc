// Draw fields of a 3D spherical detector
using namespace GeFiCa;
void draw3d()
{
   // configure detector
   Sphere3D *num = new Sphere3D;
   num->InnerR=0.5*cm;
   num->OuterR=2.5*cm;
   num->SetAverageImpurity(3e9/cm3);
   num->V0=900*volt;
   num->V1=0*volt;
   num->Dump();
   cout<<"press any key to continue"<<endl; cin.get();
   num->CalculatePotential(kSOR2); // calculate potential using SOR method

   //use analytic method from 1D sphere
   Sphere1D *ana = new Sphere1D;
   ana->InnerR=0.5*cm;
   ana->OuterR=2.5*cm;
   ana->SetAverageImpurity(3e9/cm3);
   ana->V0=900*volt;
   ana->V1=0*volt;
   ana->CalculatePotential(kAnalytic); // fill grid with analytic result

   // prepare drawing style
   gROOT->SetStyle("Plain"); // pick up a good drawing style to modify
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.02);

   // generate graphics
   TTree *tn = num->GetTree();
   tn->Draw("v:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TTree *ta = ana->GetTree();
   ta->Draw("v:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // compare numerical result to analytic calculation
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Radius [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

   TLegend *l = new TLegend(0.7,0.6,0.9,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR2","p");
   l->Draw();

   gPad->Print("s3d.png");
}
