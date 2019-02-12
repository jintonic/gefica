{
   GeFiCa::PointContactDZ *detector2 = new GeFiCa::PointContactDZ(692,505);
   detector2->Radius=3.45*GeFiCa::cm;
   detector2->Z=5.05*GeFiCa::cm;
   detector2->Rpc=0.14*GeFiCa::cm;
   detector2->Zpc=0.21*GeFiCa::cm;
   detector2->RwrapArround=3.350*GeFiCa::cm;
   detector2->TaperLength=1*GeFiCa::cm;

   TF3 *im=new TF3("f","-0.318e10+0.025e10*y");

   detector2->MaxIterations=1e5;
   detector2->Precision=1e-8;
   detector2->Csor=1.994;
   detector2->V0=500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;

   detector2->SetImpurity(im);
   
   detector2->CalculatePotential(GeFiCa::kSOR2);
   //cout<<detector2->IsDepleted()<<endl;
   detector2->SaveField("ppc2dSOR2.root");
   detector2->SaveFieldAsFieldgen("ppc2dSOR2.txt");
   //detector2->LoadField("point21dSOR23.root");

   //TCanvas * cvs=new TCanvas();
   gStyle->SetOptTitle(kFALSE);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadRightMargin(1.01);
   gStyle->SetPadLeftMargin(0.0999999999);
   gStyle->SetLabelFont(22,"XY");
   gStyle->SetLabelSize(0.06,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetTitleFont(22,"XY");
   gStyle->SetLegendFont(22);
   gStyle->SetOptStat(0);
   //gStyle->SetPalette(1);

   //hframe->Draw(); //you can set the axis att via hframe->GetXaxis()..
   TH2F* h= new TH2F("hist","",10,-3.45,3.45,10,0,5.05);
   TChain *ta = new TChain("t");
   ta->Add("ppc2dSOR2.root");
   ta->Draw("c2:c1:v>>hist","","colz");
   h->GetYaxis()->SetTitle("Thickness [cm]");
   h->GetXaxis()->SetTitle("Radius [cm]");
   //th->FillN(ta->Get("c1"),ta->Get("c2"),ta->Get("p"));
   //TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());
   //ga->GetXaxis()->SetRangeUser(-3.45,3.45);
   //ga->GetYaxis()->SetRangeUser(0,5.0);
   //ga->Draw("ap");

}
