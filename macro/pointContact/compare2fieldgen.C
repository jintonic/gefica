void compare2fieldgen()
{
   double r,z,p;

   //TH2F *hmjd = new TH2F("hmjd","",346,0,346,506,0,506);
   //ifstream ff("ev.dat");
   //for (int i = 0; i < 346; i++) {
   //  for (int j = 0; j < 506; j++) {
   //     ff>>r>>z>>p;
   //     hmjd->SetBinContent(i,j,p);
   //  } 
   //}

   TH2F *hgfc = new TH2F("hgfc","",346,0,346,506,0,506);
   ifstream fg("full.txt");
   for (int i = 0; i < 346; i++) {
     for (int j = 0; j < 506; j++) {
        fg>>r>>z>>p;
        hgfc->SetBinContent(i,j,p);
     } 
   }

   //hmjd->Draw("colz");
   //TCanvas *can = new TCanvas;
   hgfc->Draw("colz");
}

void generateField()
{
   GeFiCa::PointContactRZ *ppc = new GeFiCa::PointContactRZ(200,100);
   ppc->Radius=0.5;
   ppc->ZUpperBound=1;
   ppc->PointR=0.1;
   ppc->PointDepth=0.01;

   ppc->MaxIterations=1e6;
   ppc->Precision=1e-8;
   ppc->Csor=1.990;
   ppc->V0=0*GeFiCa::volt;
   ppc->V1=-2500*GeFiCa::volt;
   ppc->Impurity="-0.318e10+0*y";

   ppc->CalculateField(GeFiCa::kSOR2);
   ppc->SaveField("full.root");
}

void printapoint()
{
   GeFiCa::halfPointContactRZ *ppc = new GeFiCa::halfPointContactRZ;
   ppc->LoadField("half.root");
   cout<<ppc->GetPotential(0.5,0.5)<<endl;
   cout<<ppc->GetPotential(0.,0.)<<endl;
   cout<<ppc->GetPotential(0.,1)<<endl;
   cout<<ppc->GetPotential(0.,0.5)<<endl;
}

void r2t()
{
   GeFiCa::halfPointContactRZ *ppc = new GeFiCa::halfPointContactRZ(1036,506);
   ppc->LoadField("full.root");

   ofstream fo("full.txt" );
   for (int i = 0; i < 346; i++) {
     for (int j = 0; j < 506; j++) {
        fo<<1.0*i*0.5/346<<"\t"<<1.0*j*1/506<<"\t"
           <<ppc->GetPotential(1.0*i*0.5/346,1.0*j*1/506)<<endl;
      } 
   }
   fo.close(); 
}
void treedraw()
{
   TTree *t = new TTree("t","t");
   t->ReadFile("full.txt", "r:z:v");
   //t->AddFriend("t2=t","point2dSOR2.root");
   //t->Draw("z:(t2.p-v)","z!=1&r!=1&z<1&r>34.&r<34.5","");
   //TCanvas *can = new TCanvas;
   //t->Draw("r:(t2.p-v)","z>=0&z<0.2","");
   t->Draw("z:r:v","","colz");

}
