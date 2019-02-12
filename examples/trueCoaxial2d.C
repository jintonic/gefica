{
   GeFiCa::TrueCoaxial1D *detector = new GeFiCa::TrueCoaxial1D(101);
   detector->V0=0*GeFiCa::volt;
   detector->V1=2500*GeFiCa::volt;
   detector->MaxIterations=1e5;
   detector->Csor=1.95;
   detector->InnerRadius=0.14;
   detector->OuterRadius=3.45;

   GeFiCa::TrueCoaxial2D *detector2 = new GeFiCa::TrueCoaxial2D(346,505);
   detector2->InnerRadius=0.14;//-3.45;
   detector2->OuterRadius=3.45;

   TF3 *im=new TF3("f","-0.318e10+0.025e10*y");
   detector2->SetImpurity(im);


   //TF1 *im1=new TF1("f","-0.318e10+0.025e10*x",0,6.9);
   detector2->MaxIterations=1e5;
   detector2->Csor=1.994;
   detector2->V0=2500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;
   detector2->CalculatePotential(GeFiCa::kSOR2);
   detector2->SaveField("trueCoaxial2d.root");

   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("trueCoaxial1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("trueCoaxial2d.root");
   tn->Draw("v:c1*10","c2<1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("trueCoaxial1dTrue.root");
   ta->Draw("v:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   ga->SetMarkerColor(kRed);
   gn->SetMarkerStyle(6);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   ga->SetTitle(";Thickness [cm];Potential [V]");
   ga->Draw("ap");
   gn->Draw("l");
}
