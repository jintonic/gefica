{
   GeFiCa::PointContactRZ *detector2 = new GeFiCa::PointContactRZ(692,505);
   detector2->Radius=3.45;
   detector2->ZUpperBound=5.05;
   detector2->PointR=0.14;
   detector2->PointDepth=0.21;
   detector2->ContactInnerR=1.450;

   //TF2 *im=new TF2("f","-0.19175e10-0.025e10*y");
   TF3 *im=new TF3("f","-0.318e10+0.025e10*y");
   //TF1 *im1=new TF1("f","-0.318e10+0.025e10*x",0,6.9);

   detector2->MaxIterations=1e5;
   detector2->Precision=1e-8;
   detector2->Csor=1.994;
   detector2->V0=500*GeFiCa::volt;
   detector2->V1=0*GeFiCa::volt;

   //TF1 *im=new TF1("","pol1",-0.318e10,0.025e10)
   //detector2->Impurity="-0.318e10+0.025e10*y";//-0.01e10/GeFiCa::cm3);
   detector2->SetImpurity(im);
   
   detector2->CalculatePotential(GeFiCa::kSOR2);
   //cout<<detector2->IsDepleted()<<endl;
   detector2->SaveField("point2dSOR2.root");
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
   ta->Add("point2dSOR2.root");
   ta->Draw("c2:c1:p>>hist","","colz");
   h->GetYaxis()->SetTitle("Thickness [cm]");
   h->GetXaxis()->SetTitle("Radius [cm]");
   //th->FillN(ta->Get("c1"),ta->Get("c2"),ta->Get("p"));
   //TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());
   //ga->GetXaxis()->SetRangeUser(-3.45,3.45);
   //ga->GetYaxis()->SetRangeUser(0,5.0);
   //ga->Draw("ap");

/*
   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(505);
   detector->UpperBound=5.05;
   detector->V1=2500*GeFiCa::volt;
   detector->V0=0*GeFiCa::volt;
   detector->SetImpurity(-0.01e10/GeFiCa::cm3);
   detector->CalculateField(GeFiCa::kAnalytic);
   detector->SaveField("planar1dTrue.root");

  
   TCanvas * cvs=new TCanvas();
   gStyle->SetOptTitle(kTRUE);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.0999999999);
   gStyle->SetLabelFont(22,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetTitleFont(22,"XY");
   gStyle->SetLegendFont(22);
   gStyle->SetCanvasColor(kBlack);
   //cvs->SetFillColor(kBlack);

   TChain *ta = new TChain("t");
   ta->Adinputd("planar1dTrue.root");
   ta->Draw("e1:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // generate graphics
   TChain *tn = new TChain("t");
   /input/tn->Add("point2dSOR2.root");
   //tn->Draw("c2*10:p","c1<0.00&&c1>-0.05");
  // TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());
  
   TChain *ta = new TChain("t");
  // ta->Add("point2dSOR2.root");
   //ta->Draw("c2:c1:p","","colz");
   //TGraph *gn = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());
  

   
  TTree *t = new TTree("t","t");
  t->ReadFile("/home/byron/mjd_siggen/fields/p1/ev.new", "r:z:v");
  //t->AddFriend("t2=t","point2dSOR2.root");
  //t->Draw("z:(t2.p-v)","z!=1&r!=1&z<1&r>34.&r<34.5","");
  //TCanvas *can = new TCanvas;
  //t->Draw("r:(t2.p-v)","z>=0&z<0.2","");
  t->Draw("z:r:v","","colz");
  

  // TGraph *gn = new TGraph(t->GetSelectedRows(), t->GetV2(), t->GetV1());
nput
   // make final plot
   //gn->SetMarkerColor(kBlue);
   //gn->SetMarkerStyle(8);
   //gn->SetMarkerSize(0.3);
   //ga->SetLineColor(kRed);
   //gn->GetXaxis()->SetTitle("Thickness [cm]");
   //gn->GetYaxis()->SetTitle("Potential [V]");
   //gn->SetTitle("");
   
   //gn->Draw("ap");
  
   TLegend *leg = new TLegend(0.2,0.6,0.5,0.8);
   leg->SetBorderSize(0);
   //leg->AddEntry(gn,"GeFiCa","p");
   //leg->AddEntry(ga,"mjd","l");
   leg->SetTextSize(0.05);
   leg->Draw();
 //  cvs->SaveAs("pointContact2d.png");

  */ 
}
