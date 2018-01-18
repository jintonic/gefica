{
  GeFiCa::PointContactRZ *detector2 = new GeFiCa::PointContactRZ(1036,506);
  detector2->LoadField("point2dSOR2.root");

 
  std::ifstream infile("ev.new");
  std::ofstream outfile("result.txt");
  double x,y,v,anotherV;
  while (  infile>>x>>y>>v)
  {
    anotherV=detector2->GetPotential(x,y);
    outfile<<x<<"  "<<y<<"  "<<anotherV<<"  "<<v-anotherV<<endl;
  }
  infile.close();
  outfile.close();

}
