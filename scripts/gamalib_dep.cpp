#include <iostream>
#include <fstream>
#include <string>
#include <set>

/*
 * $Id: gamalib_dep.cpp,v 1.5 2003/02/22 19:40:55 cepek Exp $
 */

/*************************************************************************
 * 0.6  added directory <gnu_gama/ ... > for processing
 * 0.5  added include <iostream> for g++ 3.0.4 
 * 0.4  .o changed to .$(OBJ) for
 * 0.3  `name' not written to output
 * 0.2  added SRC make macro (2000-11-11)
 *  
 *************************************************************************/

using namespace std;

string path = "./";

void add_dep(string file, set<string>& dep)
{
  ifstream inp(file.c_str());
  if (!inp) 
    {
      file = path + file;
      inp.clear();
      inp.open(file.c_str());
      if (!inp)
        {
          cerr << "******  gamalib_dep : cannot open file " << file << endl;
          return;
        }
    }
  string line;
  while (getline(inp, line))
    {
      size_t  n, n1, n2; 
      if (string::npos == (n=line.find("#include" ))) continue;
      n1 = line.find("<gamalib/" );
      n2 = line.find("<gnu_gama/");
      if (string::npos == n1 && string::npos == n2) continue;
      if (string::npos != n1) n = n1;
      if (string::npos != n2) n = n2;
          
      string name;
      n++;
      while(line[n] != '>') name += line[n++];

      dep.insert(name);
      add_dep(name, dep);
    } 
}

int main(int argc, char* argv[])
{
  if (argc > 1)  path = string(argv[1]) + "/";

  string file, line;
  while(getline(cin, file))
    {
      set<string> dep;
      string name;
      for (string::const_iterator i=file.begin(); i!=file.end(); ++i)
        if     (*i == '/') 
          name.erase();
        else
          name += *i;

      add_dep(file, dep);
 
      
      /* # 0.3 for (string::const_iterator i=name.begin(); i!=name.end(); ++i)
       * # 0.3 {
       * # 0.3 if (*i == '.') break;
       * # 0.3 cout << *i;
       * # 0.3 }
       */
      if (file[0] == '.' && file[1] == '.' && file[2] == '/' &&
          file[3] == '.' && file[4] == '.' && file[5] == '/' &&
          file[6] == '.' && file[7] == '.' && file[8] == '/')
        {
          file.replace(0, 9, "$(SRC)");   // # 0.2
        }
      cout << ".$(OBJ) : " << file;
      for (set<string>::const_iterator i=dep.begin(); i!=dep.end(); ++i)
        cout << " $(SRC)" << *i;  // # 0.2 cout << " ../../" << *i;
      cout << endl;
    }
}











