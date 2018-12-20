{
   // draw MJD result
   TCanvas *c1 = new TCanvas;
   TTree *tm1 = new TTree("1tm","1tm");
   tm1->ReadFile("ev.new", "r:z:v:e:er:ez");
   tm1->Draw("z:r:v","","colz");
   //t->AddFriend("t2=t","point2dSOR2.root");
   //t->Draw("z:(t2.p-v)","z!=1&r!=1&z<1&r>34.&r<34.5","");
   //TCanvas *can = new TCanvas;
   //t->Draw("r:(t2.p-v)","z>=0&z<0.2","");
   //t->Draw("z:r:d","z<50","colz");
   // draw GeFiCa result
   TCanvas *c2 = new TCanvas;
   TTree *tg = new TTree("tg","tg");
   tg->ReadFile("result.txt", "r:z:v:d:e1:e2:e:de1:de2:de");
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
   
   //electricfield
}
