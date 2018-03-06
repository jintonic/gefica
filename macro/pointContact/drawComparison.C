{
   // draw MJD result
   TCanvas *c1 = new TCanvas;
   TTree *tm = new TTree("tm","tm");
   tm->ReadFile("ev.new", "r:z:v");
   tm->Draw("z:r:v","","colz");
   //t->AddFriend("t2=t","point2dSOR2.root");
   //t->Draw("z:(t2.p-v)","z!=1&r!=1&z<1&r>34.&r<34.5","");
   //TCanvas *can = new TCanvas;
   //t->Draw("r:(t2.p-v)","z>=0&z<0.2","");
   //t->Draw("z:r:d","z<50","colz");
   // draw GeFiCa result
   TCanvas *c2 = new TCanvas;
   TTree *tg = new TTree("tg","tg");
   tg->ReadFile("result.txt", "r:z:v:d");
   tg->Draw("z:r:v","","colz");

   // draw difference
   TCanvas *c3 = new TCanvas;
   tg->Draw("z:r:d","","colz");
}
