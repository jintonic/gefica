// set a very long point contact so that the field is almost the same as a
// coaxial detector. We then compare the field in the middle part of z with the
// analytic result of a true coaxial detector field
{
   // calculate potential for true coaxial 1D analyitically
   GeFiCa::TrueCoaxial1D *tc1d = new GeFiCa::TrueCoaxial1D(333);
   tc1d->V0=2500*GeFiCa::volt;
   tc1d->V1=0*GeFiCa::volt;
   tc1d->InnerRadius=0.13;
   tc1d->OuterRadius=3.45;

   tc1d->SetImpurity(-0.318e10);
   tc1d->CalculatePotential(GeFiCa::kAnalytic);
   tc1d->SaveField("tca2.root");

   // compare 
   TTree *t1 = new TTree;
   t1->ReadFile("mjd.txt", "c1:p");

   t1->Draw("p:c1");
   double *p1 = t1->GetV1();
   double *r1 = t1->GetV2();

   TChain *t2 = new TChain("t");
   t2->Add("tca2.root");

   t2->Draw("p:c1");
   double *p2 = t2->GetV1();
   double *r2 = t2->GetV2();

   const int n = t2->GetSelectedRows();
   double p[n]={0};
   for (int i=0; i<n; i++) p[i] = p1[i] - p2[i];

   TGraph *g = new TGraph(n,r1,p);
   g->Draw("ap");
}
