{
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadTopMargin(0.02);

   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D;
   detector->LowerBound=0*GeFiCa::cm;
   detector->UpperBound=1*GeFiCa::cm;
   detector->V0=0*GeFiCa::volt;
   detector->V1=2000*GeFiCa::volt;
   TF1 *im=new TF1("f","-1e10");
   detector->SetImpurity(im);
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("p0planar1dSOR2.root");

   TF1 *im=new TF1("f","1.5e10");
   detector->SetImpurity(im);
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("p1planar1dSOR2.root");

   TF1 *im=new TF1("f","3.5e10");
   detector->SetImpurity(im);
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("p3planar1dSOR2.root");
   
   TF1 *im=new TF1("f","-1.5e10");
   detector->SetImpurity(im);
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("n1planar1dSOR2.root");
   
   TF1 *im=new TF1("f","-3.5e10");
   detector->SetImpurity(im);
   detector->CalculatePotential(GeFiCa::kAnalytic);
   detector->SaveField("n3planar1dSOR2.root");

   // generate graphics
   TChain *t0 = new TChain("t");
   t0->Add("p0planar1dSOR2.root");
   t0->Draw("p:c1");
   TGraph *gr0 = new TGraph(t0->GetSelectedRows(), t0->GetV2(), t0->GetV1());

   TChain *t1 = new TChain("t");
   t1->Add("p1planar1dSOR2.root");
   t1->Draw("p:c1");
   TGraph *gp1 = new TGraph(t1->GetSelectedRows(), t1->GetV2(), t1->GetV1());

   TChain *t3 = new TChain("t");
   t3->Add("p3planar1dSOR2.root");
   t3->Draw("p:c1");
   TGraph *gp3 = new TGraph(t3->GetSelectedRows(), t3->GetV2(), t3->GetV1());

   TChain *tn1 = new TChain("t");
   tn1->Add("n1planar1dSOR2.root");
   tn1->Draw("p:c1");
   TGraph *gn1 = new TGraph(tn1->GetSelectedRows(), tn1->GetV2(), tn1->GetV1());

   TChain *tn3 = new TChain("t");
   tn3->Add("n3planar1dSOR2.root");
   tn3->Draw("p:c1");
   TGraph *gn3 = new TGraph(tn3->GetSelectedRows(), tn3->GetV2(), tn3->GetV1());

   gr0->SetTitle("");
   gr0->GetXaxis()->SetTitle("Thickness [cm]");
   gr0->GetYaxis()->SetTitle("Potential [V]");
   gr0->GetXaxis()->SetRangeUser(0,1);
   gr0->GetYaxis()->SetRangeUser(0,2000);
   gr0->GetYaxis()->SetTitleOffset(1.2);
   gr0->SetLineColor(kRed);
   gp1->SetLineColor(kBlack);
   gp3->SetLineColor(kGreen);
   gn1->SetLineColor(kBlue);
   gn3->SetLineColor(kMagenta);

   gr0->Draw("al");
   gp1->Draw("l");
   gp3->Draw("l");
   gn1->Draw("l");
   gn3->Draw("l");

   TLegend *leg = new TLegend(0.7,0.15,0.98,0.5);
   leg->SetBorderSize(0);
   leg->AddEntry(gr0,"#rho =  0/cm^{3}","l");
   leg->AddEntry(gp1,"#rho =  1.5e10/cm^{3}","l");
   leg->AddEntry(gp3,"#rho =  3.5e10/cm^{3}","l");
   leg->AddEntry(gn1,"#rho = -1.5e10/cm^{3}","l");
   leg->AddEntry(gn3,"#rho = -3.5e10/cm^{3}","l");
   leg->Draw();

   const int n=5;
   TGraph *ge[n]={0};

   t0->Draw("e1:c1");
   ge[0] = new TGraph(t0->GetSelectedRows(), t0->GetV2(), t0->GetV1());
   t1->Draw("e1:c1");
   ge[1] = new TGraph(t1->GetSelectedRows(), t1->GetV2(), t1->GetV1());
   t3->Draw("e1:c1");
   ge[2] = new TGraph(t3->GetSelectedRows(), t3->GetV2(), t3->GetV1());
   tn1->Draw("e1:c1");
   ge[3] = new TGraph(tn1->GetSelectedRows(), tn1->GetV2(), tn1->GetV1());
   tn3->Draw("e1:c1");
   ge[4] = new TGraph(tn3->GetSelectedRows(), tn3->GetV2(), tn3->GetV1());

   ge[4]->Draw("al");
   ge[4]->SetLineColor(4);
   for (int i=0; i<n-1; i++) {
      ge[i]->SetLineColor(i);
      ge[i]->Draw("l");
   }
}
