// set a very long point contact so that the field is almost the same as a
// coaxial detector. We then compare the field in the middle part of z with the
// analytic result of a true coaxial detector field
{
   // calculate potential for point contact 2D
   GeFiCa::PointContactRZ *ppc = new GeFiCa::PointContactRZ(1036,506);
   ppc->RUpperBound=3.45;
   ppc->RLowerBound=-3.45;
   ppc->ZUpperBound=5.05;
   ppc->PointBegin=-0.135;
   ppc->PointEnd=0.135;
   ppc->PointDepth=5.05;

   ppc->MaxIterations=1e6;
   ppc->Precision=1e-8;
   ppc->Csor=1.996;
   ppc->V0=2500*GeFiCa::volt;
   ppc->V1=0*GeFiCa::volt;
   ppc->Impurity="-0.318e10+0*y";
   ppc->CalculateField(GeFiCa::kSOR2);
   ppc->SaveField("pc2d.root");
   
   // calculate potential for true coaxial 1D analyitically
   GeFiCa::TrueCoaxial1D *tc1d = new GeFiCa::TrueCoaxial1D(499);
   tc1d->V0=0*GeFiCa::volt;
   tc1d->V1=2500*GeFiCa::volt;
   tc1d->InnerRadius=0.13;
   tc1d->OuterRadius=3.45;

   tc1d->Impurity="-0.318e10";
   tc1d->CalculateField(GeFiCa::kAnalytic);
   tc1d->SaveField("tca.root");

   // compare 
   TChain *t1 = new TChain("t");
   t1->Add("pc2d.root");

   t1->Draw("p:c1","c2>2.5 && c2<2.51 && c1>0.13");
   double *p1 = t->GetV1();
   double *r1 = t->GetV2();

   TChain *t2 = new TChain("t");
   t2->Add("tca.root");

   t2->Draw("p:c1");
   double *p2 = t->GetV1();
   double *r2 = t->GetV2();

   const int n = t2->GetSelectedRows();
   double p[n]={0};
   for (int i=0; i<n; i++) p[i] = p1[i] - p2[i];

   TGraph *g = new TGraph(n,r1,p);
   g->Draw("ap");
   g->GetYaxis()->SetRangeUser(-1e-9,1e-9);
}
