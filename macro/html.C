using namespace GEFICA;

void html()
{
  X x;
  THtml html;
  html.SetProductName("GEFICA");
  html.SetHomepage("htt[://github.com");
  html.SetCharset("UTF8");
  html.SetAuthorTag("Jiu Jing, Jianchen Li");
  html.SetInputDir("../core");
  html.SetOutputDir("../doc/html");
  html.MakeAll();
}

