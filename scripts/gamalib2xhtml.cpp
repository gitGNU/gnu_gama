#include <fstream>
#include <string>

/*
 * $Id: gamalib2xhtml.cpp,v 1.1 2001/12/07 11:45:44 cepek Exp $
 */	 	

int main(int argc, char* argv[])
{
  using namespace std;
  
  string libname(argv[1]);
  string file(argv[2]);
  ifstream inp(file.c_str());
  ofstream out( ( string(argv[3])+"/"+file+".html").c_str() );

  string updir;
  for (string::iterator i=file.begin(); i!=file.end(); ++i)
    if (*i == '/')
      updir += "../";
  
  out << 
    "<?xml version=\"1.0\" ?>\n"
    "<!DOCTYPE html\n"
    "     PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\"\n"
    "     \"DTD/xhtml1-strict.dtd\">\n"
    "<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
    "<head>\n"
    "<title>" << file <<  "</title>\n"
    "</head>\n"
    "<body>\n"
    "<h1><a href=\"" << updir << "classes.html\">GaMaLib</a></h1>\n"
    "<h2>" << file <<  "</h2>\n"
    "<pre>\n";

  string iline, oline, header;
  size_t n, b, e;
  while (getline(inp, iline))
    {
      oline.erase();
      for (string::iterator i=iline.begin(); i!=iline.end(); ++i)
        {
          char c = *i;
          switch (c) 
            {
            case '<': oline += "&lt;"  ; break;
            case '>': oline += "&gt;"  ; break;
            case '&': oline += "&amp;" ; break;
            default : oline += c;
            }
        }

      n = oline.find("#include");
      if(n < string::npos) 
        do 
          {
            b = oline.find("&lt;" + libname + "/");
            if (b == string::npos) break;
            e = oline.find("&gt;");
            if (e == string::npos) break;
            b += 4;
            if (e <= b) break;
            n = e - b;
            header = string(oline, b, n); 
            oline.replace(b, n, "<a href=\"" + updir + header + ".html\">"
                          + header + "</a>");
          } while (false);

      out << oline << endl;
    }


  out << "</pre>\n</body>\n</html>\n";
}






