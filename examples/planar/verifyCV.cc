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
   double d=0;
   for(int i=0;i<n-1;i++) {
      // please think about the depletion check, we cannot use fIsDepleted[i]
      if(detector->IsDepleted()&&thickness*i/(n-1)>d&&(connect[i-1])) {
         connect[i]=true;
         d=thickness*i/(n-1);
      }
   }
   std::cout<<d<<"\n";
   return d;
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
   Double_t V[n], Cn[n], Ca[n];
   for (Int_t i=1;i<n;i++) {
      V[i] = i*50;
      Cn[i] = 1/GetDepthOfDepletionNumerically(V[i],thickness);
      Ca[i] = 1/GetDepthOfDepletionAnalytically(V[i],thickness); 
   }
   TGraph *gn = new TGraph(n,V,Cn);
   TGraph *ga = new TGraph(n,V,Ca);
   gn->SetTitle(";Capacitance [Arbitrary Unit];Bias [V]");
   gn->Draw("pal");
   ga->SetMarkerStyle(kCircle);
   ga->Draw("p");
}
