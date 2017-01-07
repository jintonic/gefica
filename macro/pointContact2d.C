{/*
   GeFiCa::PointContactXY *detector2 = new GeFiCa::PointContactXY(691,505);
   detector2->XLowerBound=-3.45;
   detector2->XUpperBound=3.45;
   detector2->YUpperBound=5.05;
   detector2->PointBegin=-0.14;
   detector2->PointEnd=0.14;

   //TF2 *im=new TF2("f","-0.318e10+0.025e10*y",0,6.9,0,5.05);
   //TF1 *im1=new TF1("f","-0.318e10+0.025e10*x",0,6.9);

   detector2->MaxIterations=1e5;
   detector2->Csor=1.994;
   detector2->V0=2500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;

   //TF1 *im=new TF1("","pol1",-0.318e10,0.025e10)
   detector2->SetImpurity(-0.01e10/GeFiCa::cm3);
   detector2->CalculateField(GeFiCa::kSOR2);
   detector2->SaveField("point2dSOR2.root");*/
/*
   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(505);
   detector->UpperBound=5.05;
   detector->V1=2500*GeFiCa::volt;
   detector->V0=0*GeFiCa::volt;
   detector->SetImpurity(-0.01e10/GeFiCa::cm3);
   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("planar1dTrue.root");


   TChain *ta = new TChain("t");
   ta->Add("planar1dTrue.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

*/
   // generate graphics
   TChain *tn = new TChain("t");
   tn->Add("point2dSOR2.root");
   tn->Draw("c2:c1:p","","colz");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

/*
   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(8);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle(";Thickness [cm];Potential [V]");
   gn->Draw("ap");
   ga->Draw("l");*/
}
