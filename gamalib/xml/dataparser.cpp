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
 *  $Id: dataparser.cpp,v 1.8 2003/01/05 18:02:33 cepek Exp $
 */

#include <gamalib/xml/dataparser.h>
#include <cstring>

using namespace std;
using namespace GaMaLib;

DataParser::DataParser(std::list<DataObject*>& obs) : objects(obs)
{

  // initial parser state
  
  state = s_start;

  // implicit startElement

  for (int s=s_error; s<=s_stop; s++)
    for (int t=t_gama_data; t<=t_unknown; t++)
      {
        next[s][t] = s_error;
        func[s][t] = &DataParser::parser_error;
      }

  // implicit characterDataHandler

  for (int n=s_error; n<=s_stop; n++)
    {
      data[n] = &DataParser::white_spaces;
    }

  // implicit endElement

  for (int e=s_error; e<=s_stop; e++)
    {
      ende[e] = 0;
    }

  // .......................................................................
  next[ s_start ][ t_gama_data ]           = s_gama_data;
  func[ s_start ][ t_gama_data ]           = &DataParser::gama_data;
  // .......................................................................
  next[ s_gama_data ][ t_text ]            = s_text;
  func[ s_gama_data ][ t_text ]            = &DataParser::text;
  data[ s_text      ]                      = &DataParser::text;
  ende[ s_text      ]                      = &DataParser::text;
  // .......................................................................
  next[ s_gama_data ][ t_adj_input_data ]  = s_adj_input_data;
  func[ s_gama_data ][ t_adj_input_data ]  = &DataParser::adj_input_data;
  ende[ s_adj_input_data ]                 = &DataParser::adj_input_data;
  // .......................................................................
  next[ s_adj_input_data ][ t_sparse_mat ] = s_sparse_mat;
  func[ s_adj_input_data ][ t_sparse_mat ] = &DataParser::no_attributes;
  ende[ s_sparse_mat                     ] = &DataParser::sparse_mat;
  // .......................................................................
  next[ s_sparse_mat       ][ t_rows ]     = s_sparse_mat_rows;
  func[ s_sparse_mat       ][ t_rows ]     = &DataParser::no_attributes;
  data[ s_sparse_mat_rows  ]               = &DataParser::add_text;
  ende[ s_sparse_mat_rows  ]               = &DataParser::add_space;
  
  next[ s_sparse_mat       ][ t_cols ]     = s_sparse_mat_cols;
  func[ s_sparse_mat       ][ t_cols ]     = &DataParser::no_attributes;
  data[ s_sparse_mat_cols  ]               = &DataParser::add_text;
  ende[ s_sparse_mat_cols  ]               = &DataParser::add_space;

  next[ s_sparse_mat       ][ t_nonz ]     = s_sparse_mat_nonz;
  func[ s_sparse_mat       ][ t_nonz ]     = &DataParser::no_attributes;
  data[ s_sparse_mat_nonz  ]               = &DataParser::add_text;  
  ende[ s_sparse_mat_nonz  ]               = &DataParser::sparse_mat_nonz;

  next[ s_sparse_mat       ][ t_row  ]     = s_sparse_mat_row;
  func[ s_sparse_mat       ][ t_row  ]     = &DataParser::sparse_mat_row;
  data[ s_sparse_mat_row   ]               = &DataParser::white_spaces;

  next[ s_sparse_mat_row   ][ t_nonz ]     = s_sparse_mat_row_n;
  func[ s_sparse_mat_row   ][ t_nonz ]     = &DataParser::no_attributes;
  data[ s_sparse_mat_row_n ]               = &DataParser::add_text; 
  ende[ s_sparse_mat_row_n ]               = &DataParser::sparse_mat_row_n; 

  next[ s_sparse_mat_row   ][ t_int  ]     = s_sparse_mat_row_i;
  func[ s_sparse_mat_row   ][ t_int  ]     = &DataParser::no_attributes;
  data[ s_sparse_mat_row_i ]               = &DataParser::add_text; 
  ende[ s_sparse_mat_row_i ]               = &DataParser::add_space; 

  next[ s_sparse_mat_row   ][ t_flt  ]     = s_sparse_mat_row_f;
  func[ s_sparse_mat_row   ][ t_flt  ]     = &DataParser::no_attributes;
  data[ s_sparse_mat_row_f ]               = &DataParser::add_text; 
  ende[ s_sparse_mat_row_f ]               = &DataParser::sparse_mat_row_f; 
  // .......................................................................
}


DataParser::data_tag DataParser::tag(const char* c)
{
  switch (*c)
    {
    case 'a':
      if (!strcmp(c, "adj-input-data")) return t_adj_input_data;
      break;
    case 'c' :
      if (!strcmp(c, "cols"          )) return t_cols;
      break;
    case 'f' :
      if (!strcmp(c, "flt"           )) return t_flt;
      break;
    case 'g' :
      if (!strcmp(c, "gnu-gama-data" )) return t_gama_data;
      break;
    case 'i':
      if (!strcmp(c, "int"           )) return t_int;
      break;
    case 'n':
      if (!strcmp(c, "nonz"          )) return t_nonz;
      break;
    case 'r':
      if (!strcmp(c, "row"           )) return t_row;  // more frequent
      if (!strcmp(c, "rows"          )) return t_rows;
      break;
    case 's':
      if (!strcmp(c, "sparse-mat"    )) return t_sparse_mat;
      break;
    case 't' :
      if (!strcmp(c, "text"          )) return t_text;
      break;
    default:
      break;
    }

  error(string("### unknown tag <") + string(c) + ">");

  return t_unknown;
}


// *****************************************************************

int DataParser::parser_error(const char *cname, const char **atts)
{
  return error(string("### tag <") + string(cname) 
               + string("> cannot be used in this context"));
}

int DataParser::no_attributes(const char *cname, const char **atts)
{
  if (*atts)
    {
      return error(string("### tag <") + string(cname) 
                 + string("> cannot have any attributes"));
    }
  return 0;
}

int DataParser::white_spaces(const char* s, int len)
{
  while (len--)
    {
      if (!isspace(s[len])) return error(T_GKF_illegal_text);
    }

  return 0;
}

int DataParser::add_text(const char* s, int len)
{
  text_buffer += string(s, len);
  return 0;
}

int DataParser::add_space()
{
  text_buffer += ' ';
  return 0;
}

// ......  <gnu-gama-data>  ................................................

int DataParser::gama_data(const char *cname, const char **atts)
{
  no_attributes  (cname, atts); // will have attribute 'version' later
  return 0;
}

// ......  <text>  .........................................................

int DataParser::text(const char *cname, const char **atts)
{
  no_attributes(cname, atts);
  text_buffer.erase();
  return 0;
}

int DataParser::text(const char* s, int len)
{
  text_buffer += string(s, len);
  return 0;
}

int DataParser::text()
{
  objects.push_back( new TextDataObject(text_buffer) );
  text_buffer.erase();
  return 0;
}

// ......  <adj-input-data>  ...............................................

int DataParser::adj_input_data(const char *cname, const char **atts)
{
  no_attributes  ( cname, atts );

  tmp_sparse_mat = 0;
  tmp_block_diagonal = 0;
  tmp_vector.reset();  
  tmp_array = 0;

  return 0;
}

int DataParser::adj_input_data()
{
  AdjInputData *data = new AdjInputData;

  if (tmp_sparse_mat    ) data->set_mat(tmp_sparse_mat);
  if (tmp_block_diagonal) data->set_cov(tmp_block_diagonal);
  if (tmp_vector.dim()  ) data->set_rhs(tmp_vector);  
  if (tmp_array         ) data->set_minx(tmp_array);
  objects.push_back( new AdjInputDataObject(data) );

  return 0;
}

// ......  <sparse-mat>  ...................................................

int DataParser::sparse_mat()
{
  return 0;
}

int DataParser::sparse_mat_nonz()
{
  std::size_t  rows, cols;
  istringstream inp(text_buffer.c_str());
  if (inp >> rows >> cols >> tmp_sparse_mat_nonz)
    {
      tmp_sparse_mat = new SparseMatrix<>(tmp_sparse_mat_nonz, rows, cols);
      text_buffer.erase();
      return 0;
    }
  return error("### bad data in tags <rows> / <cols> / <nonz>");
}

int DataParser::sparse_mat_row(const char *cname, const char **atts)
{
  no_attributes(cname, atts);
  tmp_sparse_mat->new_row();
  tmp_sparse_mat_row_nonz = 0;
  return 0;
}

int DataParser::sparse_mat_row_n()
{
  istringstream inp(text_buffer.c_str());
  if (inp >> tmp_sparse_mat_row_nonz)
    {
      text_buffer.erase();
      return 0;
    }
  return error("### bad data in tag <nonz>");
}

int DataParser::sparse_mat_row_f()
{
  istringstream inp(text_buffer.c_str());
  std::size_t  indx;
  double       flt;
  if (tmp_sparse_mat_nonz-- && tmp_sparse_mat_row_nonz--)
    if (inp >> indx >> flt)
      {
        tmp_sparse_mat->add_element(flt, indx);
        text_buffer.erase();
        return 0;
      }
  return error("### bad data in tags <nonz> / <int> / <flt>");
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
    "qwerty ...\n"
    "asdfgh ...\n"
    "zxcvbn ...\n"
    "</text>\n\n"

    "<adj-input-data>\n"
    "  <sparse-mat>\n"
    "    <rows>3</rows> <cols>2</cols> <nonz>4</nonz>\n"
    "      <row> <nonz>1</nonz> <int>1</int> <flt>1.1</flt></row>\n"
    "      <row> <nonz>1</nonz> <int>2</int> <flt>2.2</flt></row>\n"
    "      <row> <nonz>2</nonz> <int>1</int> <flt>3.3</flt>      \n"
    "                           <int>2</int> <flt>4.4</flt></row>\n"
    "  </sparse-mat>\n"
    "</adj-input-data>\n"

    "<text>\n"
    "This is a DataParser demo ...\n"
    "</text>\n\n"

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

