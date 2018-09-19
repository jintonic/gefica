{
   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(101);
   detector->MaxIterations=10000;
   detector->Csor=1.95;
   detector->UpperBound=1*GeFiCa::cm;
   detector->V1=2000*GeFiCa::volt;
   detector->V0=0*GeFiCa::volt;
//   detector->SetImpurity(1e10/GeFiCa::cm3);
   detector->Impurity="0e10";
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("p0planar1dSOR2.root");

   detector->Impurity="1.5e10";
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("p1planar1dSOR2.root");

   detector->Impurity="3.5e10";
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("p3planar1dSOR2.root");
   
   detector->Impurity="-1.5e10";
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("n1planar1dSOR2.root");
   
   detector->Impurity="-3.5e10";
   detector->CalculatePotential(GeFiCa::kSOR2);
   detector->SaveField("n3planar1dSOR2.root");

   TCanvas * cvs=new TCanvas();
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetLabelFont(22,"XY");
   gStyle->SetLabelSize(0.06,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetTitleFont(22,"XY");
   gStyle->SetLegendFont(22);
   //gStyle->SetLegendTextSize(0.04);

   // generate graphics
   TChain *t0 = new TChain("t");
   t0->Add("p0planar1dSOR2.root");
   t0->Draw("p:c1");
   TGraph *g0 = new TGraph(t0->GetSelectedRows(), t0->GetV2(), t0->GetV1());

   // generate graphics
   TChain *t1 = new TChain("t");
   t1->Add("p1planar1dSOR2.root");
   t1->Draw("p:c1");
   TGraph *g1 = new TGraph(t1->GetSelectedRows(), t1->GetV2(), t1->GetV1());
   // generate graphics
   TChain *t3 = new TChain("t");
   t3->Add("p3planar1dSOR2.root");
   t3->Draw("p:c1");
   TGraph *g3 = new TGraph(t3->GetSelectedRows(), t3->GetV2(), t3->GetV1());
   // generate graphics
   TChain *tn1 = new TChain("t");
   tn1->Add("n1planar1dSOR2.root");
   tn1->Draw("p:c1");
   TGraph *gn1 = new TGraph(tn1->GetSelectedRows(), tn1->GetV2(), tn1->GetV1());
   // generate graphics
   TChain *tn3 = new TChain("t");
   tn3->Add("n3planar1dSOR2.root");
   tn3->Draw("p:c1");
   TGraph *gn3 = new TGraph(tn3->GetSelectedRows(), tn3->GetV2(), tn3->GetV1());
     // make final plot
   g1->SetMarkerColor(kBlue);
   g1->SetMarkerStyle(8);
   g1->SetMarkerSize(0.8);
   g1->SetTitle("");
   g1->GetXaxis()->SetTitle("Thickness [cm]");
   g1->GetYaxis()->SetTitle("Potential [V]");
   g1->GetYaxis()->SetTitleOffset(1.2);
   g0->SetLineColor(kRed);
   g1->SetLineColor(kBlack);
   g3->SetLineColor(kGreen);
   gn1->SetLineColor(kBlue);
   gn3->SetLineColor(kYellow);

   g0->Draw("l");
   g1->Draw("l");
   g3->Draw("l");
   gn1->Draw("l");
   gn3->Draw("l");

   TLegend *leg = new TLegend(0.2,0.7,0.5,0.8);
   leg->SetBorderSize(0);
   leg->AddEntry(g0,"rho=0","l");
   leg->AddEntry(g1,"rho=1.5e10","l");
   leg->AddEntry(g3,"rho=3.5e10","l");
   leg->AddEntry(gn1,"rho=-1.5e10","l");
   leg->AddEntry(gn3,"rho=-3.5e10","l");
   leg->SetTextSize(0.05);
   leg->Draw();
   
   //cvs->SaveAs("planar1d.png");

}
