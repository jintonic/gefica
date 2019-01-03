{
   GeFiCa::TrueCoaxial1D *detector = new GeFiCa::TrueCoaxial1D(332);
   detector->V0=0*GeFiCa::volt;
   detector->V1=2500*GeFiCa::volt;
   detector->InnerRadius=0.14;
   detector->OuterRadius=3.45;

   detector->MaxIterations=1e5;
   detector->SetImpurity(-0.318e10);
   detector->Csor=1.95;
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("trueCoaxial1d.root");
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("trueCoaxial1dTrue.root");

   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("trueCoaxial1d.root");
   tn->Draw("p:c1*10");

    TTree *ta = new TTree("t","t");
    ta->ReadFile("/home/byron/mjd_siggen/fields/p1/v.dat", "r:z:v");
    ta->Draw("v:r","z<0.1","");

   TGraph *gn = new TGraph(t->GetSelectedRows(), tn->GetV2(), tn->GetV1());


   //TChain *ta = new TChain("t");
   //ta->Add("trueCoaxial1dTrue.root");
   //ta->Draw("p:c1*10");
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
