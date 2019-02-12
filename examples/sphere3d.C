// compare numerical result to analytic calculation for 3D sphere detector
{
   GeFiCa::Sphere *sphere3d=new GeFiCa::Sphere(101,10,10);
   sphere3d->MaxIterations=1e5;
   sphere3d->Csor=1.96;
   sphere3d->V1=0*GeFiCa::volt;
   sphere3d->V0=2000*GeFiCa::volt;
   sphere3d->InnerRadius=0.5*GeFiCa::cm;
   sphere3d->OuterRadius=2.5*GeFiCa::cm;
   sphere3d->SetAverageImpurity(1e10/GeFiCa::cm3);
   sphere3d->Dump();
   cout<<"press any key to continue"<<endl;
   cin.get();

   //use analytic method from 1D sphere
   GeFiCa::Sphere1D *sphere1d=new GeFiCa::Sphere1D(101);
   sphere1d->V1=0*GeFiCa::volt;
   sphere1d->V0=2000*GeFiCa::volt;
   sphere1d->InnerRadius=0.5*GeFiCa::cm;
   sphere1d->OuterRadius=2.5*GeFiCa::cm;
   sphere1d->SetAverageImpurity(1e10/GeFiCa::cm3);

   // calculate fields with two different methods
   sphere1d->CalculatePotential(GeFiCa::kAnalytic);
   sphere1d->SaveField("sphere1dTrue.root");
   sphere3d->CalculatePotential(GeFiCa::kSOR2);
   sphere3d->SaveField("sphere3dSOR2.root");

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
   TChain *tn = new TChain("t");
   tn->Add("sphere3dSOR2.root");
   tn->Draw("v:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("sphere1dTrue.root");
   ta->Draw("v:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(kCircle);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");

   TLegend *l = new TLegend(0.7,0.6,0.9,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"SOR2","p");
   l->Draw();

   gPad->Print("sphere3d.png");
}
