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
 *  $Id: dataparser.cpp,v 1.7 2003/01/05 12:18:31 cepek Exp $
 */

#include <gamalib/xml/dataparser.h>
#include <cstring>

using namespace std;
using namespace GaMaLib;

DataParser::DataParser(std::list<DataObject*>& obs) : objects(obs)
{

  // initial parser state
  
  state = state_start;

  // implicit startElement

  for (int s=state_error; s<=state_stop; s++)
    for (int t=tag_gama_data; t<=tag_unknown; t++)
      {
        next[s][t] = state_error;
        func[s][t] = &DataParser::t_error;
      }

  // implicit characterDataHandler

  for (int n=state_error; n<=state_stop; n++)
    {
      data[n] = &DataParser::d_ws;
    }

  // endElement

  for (int e=state_error; e<=state_stop; e++)
    {
      ende[e] = 0;
    }

  // ......  gnu-gama-data  .................................
  next[ state_start ][ tag_gama_data ] = state_gama_data;
  func[ state_start ][ tag_gama_data ] = &DataParser::t_gama_data;
  // ......  text  ..........................................
  next[ state_gama_data ][ tag_text ] = state_text;
  func[ state_gama_data ][ tag_text ] = &DataParser::t_text;
  data[ state_text ] = &DataParser::d_text;
  ende[ state_text ] = &DataParser::e_text;
  // ......  adj-input-data  ................................
  next[ state_gama_data ][ tag_adj_input_data ] = state_adj_input_data;
  func[ state_gama_data ][ tag_adj_input_data ] 
    = &DataParser::t_adj_input_data;
  ende[ state_adj_input_data ] = &DataParser::e_adj_input_data;
  // ......  adj-input-data sparse-mat  .....................
  next[ state_adj_input_data ][ tag_sparse_mat ] = state_sparse_mat;
  func[ state_adj_input_data ][ tag_sparse_mat ] 
    = &DataParser::t_sparse_mat;
  ende[ state_sparse_mat ] = &DataParser::e_sparse_mat;
  //
  next[ state_sparse_mat ][ tag_sparse_mat_rows ] = state_sparse_mat_rows;
  func[ state_sparse_mat ][ tag_sparse_mat_rows ] = 
    &DataParser::t_no_attributes;
}


DataParser::data_tag DataParser::tag(const char* c)
{
  switch (*c)
    {
    case 'a':
      if (!strcmp(c, "adj-input-data")) return tag_adj_input_data;
      break;
    case 'g' :
      if (!strcmp(c, "gnu-gama-data" )) return tag_gama_data;
      break;
    case 's':
      if (!strcmp(c, "sparse-mat"    )) return tag_sparse_mat;
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

int DataParser::t_no_attributes(const char *cname, const char **atts)
{
  if (*atts)
    {
      return error(string("### tag <") + string(cname) 
                 + string("> cannot have any attributes"));
    }
  return 0;
}

int DataParser::d_ws(const char* s, int len)
{
  while (len--)
    {
      if (!isspace(s[len])) return error(T_GKF_illegal_text);
    }

  return 0;
}

// ......  gnu-gama-data  ..................................................

int DataParser:: t_gama_data(const char *cname, const char **atts)
{
  t_no_attributes  (cname, atts); // will have attribute 'version' later
  return 0;
}

// ......  Text  ...........................................................

int DataParser:: t_text(const char *cname, const char **atts)
{
  t_no_attributes  (cname, atts);
  text_buffer.erase();
  return 0;
}

int DataParser::d_text(const char* s, int len)
{
  text_buffer += string(s, len);
  return 0;
}

int DataParser::e_text()
{
  objects.push_back( new TextDataObject(text_buffer) );
  return 0;
}

// ......  AdjInputData  ...................................................

int DataParser::t_adj_input_data(const char *cname, const char **atts)
{
  t_no_attributes  ( cname, atts );

  tmp_sparse_mat = 0;
  tmp_block_diagonal = 0;
  tmp_vector.reset();  
  tmp_array = 0;

  return 0;
}

int DataParser::e_adj_input_data()
{
  AdjInputData *data = new AdjInputData;

  if (tmp_sparse_mat    ) data->set_mat(tmp_sparse_mat);
  if (tmp_block_diagonal) data->set_cov(tmp_block_diagonal);
  if (tmp_vector.dim()  ) data->set_rhs(tmp_vector);  
  if (tmp_array         ) data->set_minx(tmp_array);
  objects.push_back( new AdjInputDataObject(data) );

  return 0;
}

int DataParser::t_sparse_mat(const char *cname, const char **atts)
{
  t_no_attributes  ( cname, atts );

  text_buffer.erase();
  return 0;
}

int DataParser::e_sparse_mat()
{
  tmp_sparse_mat = new SparseMatrix<>;
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
    "</text>\n\n"

    "<adj-input-data>\n"
    "  <sparse-mat>\n"
    "  </sparse-mat>\n"
    "</adj-input-data>\n"

    "\n</gnu-gama-data>\n\n"
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
          delete *i;
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

