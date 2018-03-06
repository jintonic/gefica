{
   GeFiCa::PointContactRZ *detector2 = new GeFiCa::PointContactRZ(1036,506);
   detector2->LoadField("point2dSOR2.root");

   ifstream infile("ev.new");
   ofstream outfile("result.txt");
   double x,y,v,anotherV;
   while (infile>>x>>y>>v) {
      anotherV=detector2->GetPotential(x/10,y/10);
      outfile<<x<<"  "<<y<<"  "<<anotherV<<"  "<<v-anotherV<<endl;
   }
   infile.close();
   outfile.close();
}
