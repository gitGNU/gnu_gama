/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2004  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 * $Id: gama-g3.cpp,v 1.7 2004/01/26 19:03:09 cepek Exp $
 */

#include <fstream>
#include <iostream>
#include <gnu_gama/xml/dataparser.h>
#include <gnu_gama/list.h>
#include <gnu_gama/g3/g3_model.h>

namespace 
{
  const char* input  = 0;
  const char* output = 0;
  const char* projeq = 0;

  int error(const char* s) { std::cerr << s << "\n"; return 1; } 

  int help(int argc, char* argv[])
  {
    bool ok = argc > 1;

    for (int n=0, i=1; i<argc; i++)
      {
        const std::string a = argv[i];

        if (a == "-h" || a == "-help" || a == "--help")  
          { 
            ok = false;
            continue;
          }
        if (a == "-project-equations" || a == "--project-equations")
          {
            ++i;
            if (i < argc) 
              projeq = argv[i];
            else
              ok = false;

            continue;
          }
        
        ++n;
        if      (n == 1) input  = argv[i];
        else if (n == 2) output = argv[i];
        else 
          {
            ok = false;
          }
      }

    if (ok) return 0;

    std::cerr << 
      "\n"
      "Usage:  gama-g3  [ options ] input  [ output ] \n\n"

      " input      xml data file name\n"
      " output     optional output data file name\n\n"

      " -project-equations file"
      "     optional output of project equations in XML\n"

      "\n"
      " -h         help (this text)\n"

      "\n";

    return 1;
  }
  
  GNU_gama::g3::Model* get_xml_input(const char* file)
  {
    using namespace GNU_gama::g3;

    std::ifstream input(file);
    if (!input) return 0;

    GNU_gama::List<GNU_gama::DataObject::Base*> objects;
    GNU_gama::DataParser parser(objects);

    try 
      {
        std::string text;
        while (std::getline(input, text))
          {
            parser.xml_parse(text.c_str(), text.length(), 0);
            parser.xml_parse("\n", 1, 0);
          }
        parser.xml_parse("", 0, 1);
      }
    catch(const GNU_gama::Exception::parser& p)
      {
        std::cerr << "\nXML parser error on line " << p.line 
                  << " of input data  "
                  << "\t(error code " << p.error_code << ")\n"
                  << p.str << "\n\n";
        return 0;
      }
    catch(...)
      {
        error("catch ... ");
        return 0;
      }

    GNU_gama::g3::Model* model = 0;
    GNU_gama::List<GNU_gama::DataObject::Base*>::iterator i = objects.begin();
    GNU_gama::List<GNU_gama::DataObject::Base*>::iterator e = objects.end();
    while (i != e)
      {
         if (GNU_gama::DataObject::g3_model* m
            = dynamic_cast<GNU_gama::DataObject::g3_model*>(*i))
          {
            if (model) delete model;   // this should never happen
            model =  m->model;
          }

         // std::cerr << (*i)->xml() << std::endl;
         delete *i;
         ++i;
      }

    return model;
  }
}

// ----------------------------------------------------------------------

int main(int argc, char* argv[])
{
  if (help(argc, argv)) return 1;

  using namespace GNU_gama::g3;

  Model* model = get_xml_input(input);
  if (model == 0) return 0; // ??? error("error on reading XML input data");
  
  if (model) 
    {
      using namespace std;
      cerr.precision(12);
      Model::ObservationData::iterator i = model->obsdata.begin();
      Model::ObservationData::iterator e = model->obsdata.end();
      while (i != e)
        {
          cerr << "* ";
          if (Distance *d = dynamic_cast<Distance*>(*i))
            {
              cerr << " distance : from = "
                   << d->from
                   << "  to = "
                   << d->to
                   << "  val = "
                   << d->obs();
              
            }
          cerr << "\n";
          ++i;
        }

      model->update_linearization();

      if (projeq) 
        {
          std::ofstream out(projeq);
          out.precision(16);
          out << GNU_gama::DataParser::xml_start;
          model->write_xml_adjustment_input_data(out);
          out << GNU_gama::DataParser::xml_end;
        }
    }
  
  delete model;

  return 0;
}
