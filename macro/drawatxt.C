{

  TTree *t = new TTree("t","t");
  t->ReadFile("result.txt", "r:z:v:d");
  //t->AddFriend("t2=t","point2dSOR2.root");
  //t->Draw("z:(t2.p-v)","z!=1&r!=1&z<1&r>34.&r<34.5","");
  //TCanvas *can = new TCanvas;
  //t->Draw("r:(t2.p-v)","z>=0&z<0.2","");
  t->Draw("z:r:d","z<35","colz");




}
