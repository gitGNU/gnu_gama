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
 *  $Id: dataparser.cpp,v 1.9 2003/01/06 17:44:12 cepek Exp $
 */

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
    "    <rows>3</rows> <cols>2</cols> <nonz>999</nonz>\n"
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

#else
// #########################################################################


#include <gamalib/xml/dataparser.h>
#include <cstring>

using namespace std;
using namespace GaMaLib;

DataParser::DataParser(std::list<DataObject*>& obs) : objects(obs)
{
  int s, t, n, z;

  // initial parser state and implicit handlers
  
  state = s_start;

  for (s=s_error; s<=s_stop; s++)
    {
      for (t=0; t<=t_unknown; t++)
        {
          next[s][t] = s_error;
          stag[s][t] = &DataParser::parser_error;
        }
      after[s] = s_error;
      data [s] = &DataParser::white_spaces;
      etag [s] = &DataParser::end_tag;
    }


  // .....  <gnu-gama-data>  .........................................

  s = s_start;                    
  t = t_gama_data;        
  n = s_gama_data;

  next [s][t] = n;                  stag[s][t] = &DataParser::gama_data;  
  after[n]    = s_stop;

  // .....  <text>  ..................................................
 
  s = s_gama_data;                
  t = t_text;
  n = s_text;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_gama_data;        data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::text;
 
  // .....  <adj-input-data>  ........................................

  s = s_gama_data;
  t = t_adj_input_data;
  n = s_adj_input_data;

  next [s][t] = n;                  stag[s][t] = &DataParser::adj_input_data;
  after[n]    = s_gama_data;        etag[n]    = &DataParser::adj_input_data;

  // .....  <sparse-mat>  ............................................

  s = s_adj_input_data;
  t = t_sparse_mat;
  n = s_sparse_mat_1;
  z = s_sparse_mat_4;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[z]    = s_adj_input_data;   etag[z]    = &DataParser::sparse_mat;

  s = s_sparse_mat_1;
  t = t_rows;
  n = s_sparse_mat_rows;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_sparse_mat_2;     data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::append_sp;

  s = s_sparse_mat_2;
  t = t_cols;
  n = s_sparse_mat_cols;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_sparse_mat_3;     data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::append_sp; 

  s = s_sparse_mat_3;
  t = t_nonz;
  n = s_sparse_mat_nonz;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_sparse_mat_4;     data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::sparse_mat_nonz;

  s = s_sparse_mat_4;
  t = t_row;
  n = s_sparse_mat_row_1;
  z = s_sparse_mat_row_2;

  next [s][t] = n;                  stag[s][t] = &DataParser::sparse_mat_row;
  after[z]    = s_sparse_mat_4;     etag[z]    = &DataParser::sparse_mat_row;

  s = s_sparse_mat_row_1;
  t = t_nonz;
  n = s_sparse_mat_row_n;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_sparse_mat_row_2; data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::sparse_mat_row_n;

  s = s_sparse_mat_row_2;
  t = t_int;
  n = s_sparse_mat_row_i;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_sparse_mat_row_3; data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::append_sp; 

  s = s_sparse_mat_row_3;
  t = t_flt;
  n = s_sparse_mat_row_f;

  next [s][t] = n;                  stag[s][t] = &DataParser::start_tag;
  after[n]    = s_sparse_mat_row_2; data[n]    = &DataParser::add_text;
                                    etag[n]    = &DataParser::sparse_mat_row_f;

  // .................................................................
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

int DataParser::parser_error(const char *name, const char **atts)
{
  return error(string("### tag <") + string(name) 
               + string("> cannot be used in this context"));
}


int DataParser::start_tag(const char *name, const char **atts)
{
  no_attributes(name, atts);
  return (state = next[state][tag(name)]);
}

int DataParser::end_tag(const char *name)
{
  return (state = after[state]);
}

int DataParser::no_attributes(const char *name, const char **atts)
{
  if (*atts)
    {
      return error(string("### tag <") + string(name) 
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

int DataParser::append_sp(const char *name)
{
  text_buffer += ' ';
  return end_tag(name);
}

// ......  <gnu-gama-data>  ................................................

int DataParser::gama_data(const char *name, const char **atts)
{
  no_attributes  (name, atts);   // will have attribute 'version' later?
  state = next[state][tag(name)];
  return 0;
}

// ......  <text>  .........................................................

int DataParser::text(const char* name)
{
  objects.push_back( new TextDataObject(text_buffer) );
  text_buffer.erase();
  return end_tag(name);
}

// ......  <adj-input-data>  ...............................................

int DataParser::adj_input_data(const char *name, const char **atts)
{
  no_attributes( name, atts );
  state = next[state][tag(name)];

  adj_sparse_mat = 0;
  adj_block_diagonal = 0;
  adj_vector.reset();  
  adj_array = 0;

  return 0;
}

int DataParser::adj_input_data(const char *name)
{
  AdjInputData *data = new AdjInputData;

  if (adj_sparse_mat    ) data->set_mat(adj_sparse_mat);
  if (adj_block_diagonal) data->set_cov(adj_block_diagonal);
  if (adj_vector.dim()  ) data->set_rhs(adj_vector);  
  if (adj_array         ) data->set_minx(adj_array);
  objects.push_back( new AdjInputDataObject(data) );

  return end_tag(name);;
}

// ......  <sparse-mat>  ...................................................

int DataParser::sparse_mat(const char *name)
{
  return end_tag(name);
}

int DataParser::sparse_mat_nonz(const char *name)
{
  std::size_t  rows, cols;
  istringstream inp(text_buffer.c_str());
  if (inp >> rows >> cols >> adj_sparse_mat_nonz)
    {
      adj_sparse_mat = new SparseMatrix<>(adj_sparse_mat_nonz, rows, cols);
      text_buffer.erase();
      return end_tag(name);
    }
  return error("### bad data in tags <rows> / <cols> / <nonz>");
}

int DataParser::sparse_mat_row(const char *name, const char **atts)
{
  no_attributes(name, atts);
  state = next[state][tag(name)];

  adj_sparse_mat->new_row();
  adj_sparse_mat_row_nonz = 0;
  return 0;
}

int DataParser::sparse_mat_row(const char *name)
{
  return end_tag(name);
}

int DataParser::sparse_mat_row_n(const char *name)
{
  istringstream inp(text_buffer.c_str());
  if (inp >> adj_sparse_mat_row_nonz)
    {
      text_buffer.erase();
      return end_tag(name);
    }
  return error("### bad data in tag <nonz>");
}

int DataParser::sparse_mat_row_f(const char *name)
{
  istringstream inp(text_buffer.c_str());
  std::size_t  indx;
  double       flt;
  if (adj_sparse_mat_nonz-- && adj_sparse_mat_row_nonz--)
    if (inp >> indx >> flt)
      {
        adj_sparse_mat->add_element(flt, indx);
        text_buffer.erase();
        return end_tag(name);
      }
  return error("### bad data in tags <nonz> / <int> / <flt>");
}


#endif
