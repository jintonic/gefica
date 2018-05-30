{
   // calculate fields
   GeFiCa::Planar1D *detector = new GeFiCa::Planar1D(11);
   detector->MaxIterations=10;
   detector->Csor=1.95;
   detector->UpperBound=0.175*GeFiCa::cm;
   detector->V1=1000*GeFiCa::volt;
   detector->V0=0*GeFiCa::volt;
   //detector->SetImpurity(1e10/GeFiCa::cm3);
//   detector->Impurity="1e10";

   //detector->CalculatePotential(GeFiCa::kSOR2);
   detector->CalculatePotential(GeFiCa::kCG);
   detector->SaveField("planar1dSOR2.root");
  /* 
     
 //  detector->Impurity="0";
   for (int i=0;i<10;i++)
   {
      detector->MaxIterations=30*i;
      detector->CalculatePotential(GeFiCa::kSOR2);
      detector->SaveField("planar1dSOR2.root");
      TCanvas * cvs=new TCanvas();
      TChain *tn = new TChain("t");
 
      tn->Add("planar1dSOR2.root");
      tn->Draw("p:c1");

  
   }
  



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
   TChain *tn = new TChain("t");
   tn->Add("planar1dSOR2.root");
   tn->Draw("p:c1");
   TGraph *gn = new TGraph(tn->GetSelectedRows(), tn->GetV2(), tn->GetV1());

   TChain *ta = new TChain("t");
   ta->Add("planar1dTrue.root");
   ta->Draw("p:c1");
   TGraph *ga = new TGraph(ta->GetSelectedRows(), ta->GetV2(), ta->GetV1());

   // make final plot
   gn->SetMarkerColor(kBlue);
   gn->SetMarkerStyle(8);
   gn->SetMarkerSize(0.8);
   ga->SetLineColor(kRed);
   gn->SetTitle("");
   gn->GetXaxis()->SetTitle("Thickness [cm]");
   gn->GetYaxis()->SetTitle("Potential [V]");
   gn->GetYaxis()->SetTitleOffset(1.2);

   gn->Draw("ap");
   ga->Draw("l");

   TLegend *leg = new TLegend(0.2,0.6,0.5,0.8);
   leg->SetBorderSize(0);
   leg->AddEntry(gn,"SOR2","p");
   leg->AddEntry(ga,"Analyic","l");
   leg->SetTextSize(0.05);
   leg->Draw();
   
   //cvs->SaveAs("planar1d.png");
*/
}
