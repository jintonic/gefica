using namespace GeFiCa;
/**
 * \file compare2RhoZ.cc
 * \example compare2RhoZ.cc
 * \brief Compare PointConactDZ to PointContactRhoZ
 */
void compare2RhoZ()
{
   // calculate potential for point contact 2D
   GeFiCa::PointContactRZ *ppc = new GeFiCa::PointContactRZ(1036,506);
   ppc->Radius=3.45;
   ppc->ZUpperBound=5.05;
   ppc->PointR=0.135;
   ppc->PointDepth=1.05;

   ppc->MaxIterations=1e6;
   ppc->Precision=1e-8;
   ppc->Csor=1.996;
   ppc->V0=2500*GeFiCa::volt;
   ppc->V1=0*GeFiCa::volt;
   ppc->SetImpurity(-0.318e10);
   ppc->CalculatePotential(GeFiCa::kSOR2);
   //ppc->SaveField("pc2d.root");
   
   
   //calculate for half of point contact 2D
   GeFiCa::halfPointContactRZ *hppc = new GeFiCa::halfPointContactRZ(518,506);
   hppc->Radius=3.45;
   hppc->ZUpperBound=5.05;
   hppc->PointR=0.135;
   hppc->PointDepth=1.05;

   hppc->MaxIterations=1e6;
   hppc->Precision=1e-8;
   hppc->Csor=1.993;
   hppc->V0=2500*GeFiCa::volt;
   hppc->V1=0*GeFiCa::volt;
   hppc->SetImpurity(-0.318e10);
   hppc->CalculatePotential(GeFiCa::kSOR2);
   //hppc->SaveField("hpc2d.root");
   
     // compare 
   TChain *t1 = new TChain("t");
   t1->Add("pc2d.root");
 //  t1->AddFriend("t2=t","hpc2d.root");
 //  t1->Draw("c1:c2:p-t2.p","c1>0","colz");


   t1->Draw("p:c1","c2>0.1 && c2<0.11 && c1>0");
   double *p1 = t1->GetV1();
   double *r1 = t1->GetV2();

   TChain *t2 = new TChain("t");
   t2->Add("hpc2d.root");

   t2->Draw("p:c1","c2>0.1 && c2<0.11 && c1>0");
   double *p2 = t2->GetV1();
   double *r2 = t2->GetV2();


    const int n = t1->GetSelectedRows();
   double p[n]={0};
   for (int i=0; i<n; i++) {
      cout<<setprecision(20)<<"p1["<<i<<"]="<<p1[i]<<", p2["<<i<<"]="<<p2[i]<<endl;
      p[i] = p1[i] - p2[i];
   }

   TGraph *g = new TGraph(n,r1,p);
   g->Draw("ap");
   //g->GetYaxis()->SetRangeUser(-1e-9,1e-9);
}
