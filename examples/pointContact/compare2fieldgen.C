using namespace GeFiCa;
/**
 * \file compare2fieldgen.cc
 * \example pointContact/compare2fieldgen.cc
 * \brief Compare PointContactDZ to fieldgen
 */
void drawresult()
{
   // draw MJD result
   TCanvas *c1 = new TCanvas;
   TTree *tm1 = new TTree("tm1","tm1");
   tm1->ReadFile("ev.dat", "r:z:v:e:er:ez");
   tm1->Draw("z:r:v","","colz");
   //t->AddFriend("t2=t","ppcSOR2.root");
   //t->Draw("z:(t2.p-v)","z!=1&r!=1&z<1&r>34.&r<34.5","");
   //TCanvas *can = new TCanvas;
   //t->Draw("r:(t2.p-v)","z>=0&z<0.2","");
   //t->Draw("z:r:d","z<50","colz");
   // draw GeFiCa result
   TCanvas *c2 = new TCanvas;
   TTree *tg = new TTree("tg","tg");
   tg->ReadFile("gVSf.txt", "r:z:v:d:e1:e2:e:de1:de2:de");
   tg->Draw("z:r:v","","colz");
   TCanvas *c4 = new TCanvas;
   
   c4->SetFillColor(kBlue);
   
   tg->Draw("z:r:e","","colz");
   TCanvas *c5 = new TCanvas;
   tg->Draw("z:r:e1","","colz");
   TCanvas *c6 = new TCanvas;
   tg->Draw("z:r:e2","","colz");

   // draw difference
   TCanvas *c3 = new TCanvas;
   tg->Draw("z:r:d","","colz");
   TCanvas *c7 = new TCanvas;
   tg->Draw("z:r:de1","","colz");
   TCanvas *c8 = new TCanvas;
   tg->Draw("z:r:de2","","colz");
   TCanvas *c9 = new TCanvas;
   
   c9->SetFillColor(kBlue);
   tg->Draw("z:r:de","","colz");
}

//______________________________________________________________________________
//
void compare2fieldgen(const char *gefica="ppc2dSOR2.root",
      const char *fieldgen="ev.dat")
{
   TFile *inrootfile=new TFile("rcpc.root","READ");
   GeFiCa::PointContactDZ *detector2;
   infile ->GetObject("rcpc",inrootfile);

   ifstream infile("ev.dat"); if (!infile.is_open()) exit(-1);
   ofstream outfile("gVSf.txt");

   double x,y,v,er,ez,e,anotherV,E1,E2,E,sizeofr,sizeofz;
   while (infile>>x>>y>>v>>e>>er>>ez) {
      sizeofr=x;
      sizeofz=y;
      anotherV=detector2->GetPotential(x/10,y/10);
      E1=detector2->GetE1(x/10,y/10);
      E2=detector2->GetE2(x/10,y/10);
      E=sqrt(E1*E1+E2*E2);
      outfile<<x<<"  "<<y<<"  "<<anotherV<<"  "<<v-anotherV<<"  "<<E1<<"  "<<E2<<"  "<<E<<"  "<<er-E1<<"  "<<ez-E2<<"  "<<e-E<<endl;
   }
   infile.close();
   outfile.close();

   drawresult();
}
