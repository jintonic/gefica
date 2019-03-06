using namespace GeFiCa;
// Compare numerically calculated C-V curve with analytic calculation
// Since C=epsilon*A/d, C-V curve is basically 1/d-V curve
double GetDepthOfDepletionNumerically(double voltage, double thickness)
{
   const int n=2001;
   // calculate fields
   Planar1D *detector = new Planar1D(n);
   detector->Thickness=thickness*cm;
   detector->V1=0*volt;
   detector->V0=voltage*volt;

   TF3 *im=new TF3("f","-1e10");
   detector->SetImpurity(im);
   bool *connect=new bool[n];
   for(int i=0;i<n;i++) {
      connect[i]=false;
   }
   connect[0]=true;

   detector->CalculatePotential(kSOR2);
   
   //search for first V0
   double begin=0;
   double end=thickness;
   double mid=(begin+end)/2;
   while(begin<end)
   {
      mid=(begin+end)/2;
      if(detector->GetV(mid)>0)
      {
         begin=mid+1e-5;
      }
      if(detector->GetV(mid)<=0)
      {
         end=mid-1e-5;
      }
   }

   std::cout<<mid<<"\n";
   return mid;
}
double GetDepthOfDepletionbuildinNumerically(double voltage, double thickness)
{
   const int n=2001;
   // calculate fields
   Planar1D *detector = new Planar1D(n);
   detector->Thickness=thickness*cm;
   detector->V1=0*volt;
   detector->V0=voltage*volt;

   TF3 *im=new TF3("f","-1e10");
   detector->SetImpurity(im);
   detector->CalculatePotential(kSOR2);
   
   return detector->GetC();
}
//______________________________________________________________________________
//
double GetDepthOfDepletionAnalytically(double voltage, double thickness)
{
   double impurity=-1e10*Qe/cm3;
   double d=TMath::Sqrt(-2*epsilon*voltage/impurity);
   if (d>thickness)return thickness;
   //voltage=ax^2+c2x+c1
   //c1=0, when voltage=0 at x=0
   //c2=-ad, when E=dV/dx=0 at x=0, where just depleted
   //a is impurity/epsilon
   //voltage=-ax^2/2, solve voltage when x=d
   return d;
}
//______________________________________________________________________________
//
void verifyCV()
{
   // please don't draw diff, instead, overlay two curves on top of each other
   double thickness=1*cm;
   const int n=15;
   Double_t V[n], Cn[n], Ca[n],Cnb[n];
   for (Int_t i=0;i<n;i++) {
      V[i] = (i+1)*50;
      Cn[i] = 1/GetDepthOfDepletionNumerically(V[i],thickness);
      Cnb[i] = GetDepthOfDepletionbuildinNumerically(V[i],thickness)/pF;
      Ca[i] = 1/GetDepthOfDepletionAnalytically(V[i],thickness); 
   }
   for(int i=0;i<n;i++)
   {
      Cnb[i]=Cnb[i]/Cnb[n-1];
   }
   // prepare drawing style
   gROOT->SetStyle("Plain"); // pick up a good drawing style to modify
   gStyle->SetLegendBorderSize(0);
   gStyle->SetLegendFont(132);
   gStyle->SetLabelFont(132,"XY");
   gStyle->SetTitleFont(132,"XY");
   gStyle->SetLabelSize(0.05,"XY");
   gStyle->SetTitleSize(0.05,"XY");
   gStyle->SetPadRightMargin(0.01);
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadTopMargin(0.02);
   TGraph *gn = new TGraph(n,V,Cn);
   TGraph *gnb = new TGraph(n,V,Cnb);
   TGraph *ga = new TGraph(n,V,Ca);
   gn->SetTitle(";Bias [V];Capacitance [Arbitrary Unit]");
   gnb->SetMarkerStyle(kCircle);
   ga->SetMarkerStyle(7);
   //gnb->SetMarkerStyle(kCircle);
   gnb->Draw("ap");
   gn->Draw("*");
   ga->Draw("l");
   TLegend *l = new TLegend(0.5,0.6,0.8,0.8);
   l->AddEntry(ga,"Analytic","l");
   l->AddEntry(gn,"numerical","*");
   l->AddEntry(gnb,"buildin numerical","p");
   l->Draw();
}
