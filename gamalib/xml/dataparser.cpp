/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU GaMa / GaMaLib C++ Library.
    
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
 *  $Id: dataparser.cpp,v 1.5 2003/01/03 17:54:06 cepek Exp $
 */

#include <gamalib/xml/dataparser.h>
#include <cstring>

using namespace std;
using namespace GaMaLib;

DataParser::DataParser(std::list<DataObject*>& obs) : objects(obs)
{
  // startElement

  for (int s=state_error; s<=state_stop; s++)
    for (int t=tag_gama_data; t<=tag_unknown; t++)
      {
        fun[s][t] = &DataParser::t_error;
      }

  fun[state_start    ][tag_gama_data] = &DataParser::t_gama_data;
  fun[state_gama_data][tag_text     ] = &DataParser::t_text;


  // endElement

  next[state_error    ] = state_error; 
  next[state_start    ] = state_error;
  next[state_gama_data] = state_stop;
  next[state_text     ] = state_gama_data;
  next[state_stop     ] = state_error;  


  // characterDataHandler

  for (int n=state_error; n<=state_stop; n++)
    {
      dat[n] = &DataParser::d_ws;
    }

  dat[state_text] = &DataParser::d_text;


  // initial parser state
  
  state = state_start;
}




DataParser::data_tag DataParser::tag(const char* c)
{
  switch (*c)
    {
    case 'g' :
      if (!strcmp(c, "gnu-gama-data" )) return tag_gama_data;
      break;
    case 't' :
      if (!strcmp(c, "text"          )) return tag_text;
      break;
    default:
      break;
    }

  error(string("### unknown tag <") + string(c) + ">");

  return tag_unknown;
}


// *****************************************************************

int DataParser::t_error(const char *cname, const char **atts)
{
  return error(string("### tag <") + string(cname) 
               + string("> cannot be used in this context"));
}

int DataParser:: t_gama_data(const char *cname, const char **atts)
{
  state = state_gama_data;
  return 0;
}

int DataParser:: t_text(const char *cname, const char **atts)
{
  state = state_text;
  objects.push_back( new TextDataObject );
  return 0;
}


// *****************************************************************

int DataParser::d_ws(const char* s, int len)
{
  while (len--)
    {
      if (!isspace(s[len])) return error(T_GKF_illegal_text);
    }

  return 0;
}

int DataParser::d_text(const char* s, int len)
{
  if (TextDataObject* tdo = dynamic_cast<TextDataObject*>(objects.back()) )
    tdo->text += string(s, len);

  return 0;
}


// #########################################################################

#ifdef GaMaLib_DataParser_demo

#include <gamalib/xml/dataparser.h>
#include <cstring>
#include <iostream>


int main()
{
  using namespace std;
  using namespace GaMaLib;
  
  const char* xml_input_data = 

    "<?xml version=\"1.0\" ?>\n"
    "<!DOCTYPE gnu-gama-data SYSTEM \"gnu-gama-data.dtd\">\n\n"

    "<gnu-gama-data>\n\n"
    "<text>\n"
    "This is a DataParser demo ...\n"
    "</text>\n\n"

    "<text>\n"
    "qwerty. ..\n"
    "asdfgh ...\n"
    "zxcvbn ...\n"
    "</text>\n"
    "\n</gnu-gama-data>\n"
    ;

  try 
    {
      set_gama_language(en);

      list<DataObject*> objects;
      DataParser dp(objects);
      dp.xml_parse(xml_input_data, strlen(xml_input_data), 1);
      
      cout << DataObject::xml_begin();

      for (list<DataObject*>::const_iterator i=objects.begin(); 
           i!=objects.end(); ++i)
        {
          cout << (*i)->xml();
        }

      cout << DataObject::xml_end();
    }
  catch(ParserException e)
    {
      cout << "\nOn line " << e.line 
           << " --- exception : <"   << e.text 
           << ">  --- error code : " << e.error_code << endl;
      return 1;
    }
  catch(Exception g)
    {
      cout << "\nGaMaLib Exception : " << g.text << endl;
      return 2;
    }
  catch(...)
    {
      cout << "\nUnknown exception\n";
      return 3;
    }

  return 0;
}


#endif

