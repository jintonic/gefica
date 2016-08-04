using namespace GEFICA;

void html()
{
  X x;
  THtml html;
  html.SetProductName("GEFICA");
  html.SetHomepage("http://github.com/jintonic/gefica");
  html.SetCharset("UTF-8");
  html.SetAuthorTag("Jing Liu, Jianchen Li");
  html.SetInputDir("../core");
  html.SetOutputDir("../doc/html");
  html.MakeAll();
}

