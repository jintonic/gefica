{
  GeFiCa::Sphere1D *detector=new GeFiCa::Sphere1D(101);
  detector->MaxIterations=1e5;
  detector->Csor=1.9;
  detector->V1=0*GeFiCa::volt;
  detector->V0=2000*GeFiCa::volt;
  detector->SetImpurity(1e10/GeFiCa::cm3);
  detector->CalculateField(GeFiCa::kAnalytic);
  detector->SaveField("sphere1dTrue.root");

  GeFiCa::Sphere *detector2=new GeFiCa::Sphere(101,10,10);
  detector2->MaxIterations=1e5;
  detector2->Csor=1.999;
  detector2->V1=0*GeFiCa::volt;
  detector2->V0=2000*GeFiCa::volt;
  //detector2->SetImpurity(1e10/GeFiCa::cm3);
  detector2->CalculateField(GeFiCa::kSOR2);
  detector2->SaveField("sphere.root");


   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("sphere.root");
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
