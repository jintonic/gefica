{
   GeFiCa::TrueCoaxial1D *detector = new GeFiCa::TrueCoaxial1D(101);
   detector->V1=0*GeFiCa::volt;
   detector->V0=3000*GeFiCa::volt;
   detector->MaxIterations=1e5;
   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->Csor=1.95;


 GeFiCa::TrueCoaxial2D2 *detector2 = new GeFiCa::TrueCoaxial2D2(691,505);

  detector2->XLowerBound=-3.45;
   detector2->OuterRadius=3.45;
   detector2->=5.05;

   //TF2 *im=new TF2("f","-0.318e10+0.025e10*y",0,6.9,0,5.05);
   //TF1 *im1=new TF1("f","-0.318e10+0.025e10*x",0,6.9);

   detector2->MaxIterations=1e5;
   detector2->Csor=1.994;
   detector2->V0=2500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;

   //TF1 *im=new TF1("","pol1",-0.318e10,0.025e10)
   detector2->SetImpurity(-0.01e10/GeFiCa::cm3);
   detector2->CalculateField(GeFiCa::kSOR2);
   detector2->SaveField("trueCoaxial2d.root");

   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("trueCoaxial1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("trueCoaxial2d.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("trueCoaxial1dTrue.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(6);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   ga->Draw("al");
   gn->Draw("p");

}
