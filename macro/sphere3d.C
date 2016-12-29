{
   GeFiCa::Sphere1D *sphere1d=new GeFiCa::Sphere1D(101);
   sphere1d->MaxIterations=1e5;
   sphere1d->Csor=1.9;
   sphere1d->V1=0*GeFiCa::volt;
   sphere1d->V0=2000*GeFiCa::volt;
   sphere1d->InnerRadius=0.5*GeFiCa::cm;
   sphere1d->OuterRadius=2.5*GeFiCa::cm;
   sphere1d->SetImpurity(1e10/GeFiCa::cm3);
   sphere1d->CalculateField(GeFiCa::kAnalytic);
   sphere1d->SaveField("sphere1dTrue.root");

   GeFiCa::Sphere *sphere3d=new GeFiCa::Sphere(101,10,10);
   sphere3d->MaxIterations=1e5;
   sphere3d->Csor=1.999;
   sphere3d->V1=0*GeFiCa::volt;
   sphere3d->V0=2000*GeFiCa::volt;
   sphere3d->InnerRadius=0.5*GeFiCa::cm;
   sphere3d->OuterRadius=2.5*GeFiCa::cm;
   sphere3d->SetImpurity(1e10/GeFiCa::cm3);
   sphere3d->CalculateField(GeFiCa::kSOR2);
   sphere3d->SaveField("sphere3dSOR2.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("sphere3dSOR2.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("sphere1dTrue.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(6);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");
}
