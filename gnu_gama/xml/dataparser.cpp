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
 *  $Id: dataparser.cpp,v 1.5 2003/05/17 17:07:08 cepek Exp $
 */

// #########################################################################
#ifdef GNU_Gama_DataParser_demo

#include <gnu_gama/xml/dataparser.h>
#include <cstring>
#include <iostream>


int main()
{
  using namespace std;
  using namespace GNU_gama;
  
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
    "  </sparse-mat>\n\n"

    "  <block-diagonal>\n"
    "    <blocks>1</blocks> <nonz>678</nonz>\n"
    "      <block> <dim>3</dim> <width>0</width>\n"
    "      <flt>1.01</flt>"
    "      <flt>2.01</flt>"
    "      <flt>3.01</flt>"
    "      </block>\n"
    "  </block-diagonal>\n\n"

    "     <vector>\n"
    "       <dim>3</dim>\n"
    "         <flt>10.1</flt>\n"
    "         <flt>10.2</flt>\n"
    "         <flt>10.3</flt>\n"
    "     </vector>\n\n"

    "     <array>\n"
    "       <dim>2</dim>\n"
    "         <int>1</int>\n"
    "         <int>2</int>\n"
    "     </array>\n"

    "</adj-input-data>\n"

    "<text>\n"
    "This is a DataParser demo ...\n"
    "</text>\n\n"
    
    "\n</gnu-gama-data>\n\n"
    ;

  try 
    {
      set_gama_language(en);

      list<DataObject::Base*> objects;
      DataParser dp(objects);
      dp.xml_parse(xml_input_data, strlen(xml_input_data), 1);
      
      cout << DataObject::xml_begin();

      for (list<DataObject::Base*>::const_iterator i=objects.begin(); 
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


#include <gnu_gama/xml/dataparser.h>
#include <gnu_gama/g3/g3_cluster.h>
#include <cstring>

using namespace std;
using namespace GNU_gama;

DataParser::~DataParser()
{
  delete mg3;

  delete adj_sparse_mat;
  delete adj_block_diagonal;
  delete adj_array;
}


DataParser::DataParser(List<DataObject::Base*>& obs) : objects(obs)
{
  mg3 = 0;

  adj_sparse_mat = 0;
  adj_block_diagonal = 0;
  adj_array = 0;

  // initial parser state and implicit handlers
  
  state = s_start;

  for (int s=s_error; s<=s_stop; s++)
    {
      for (int t=0; t<=t_unknown; t++)
        {
          next[s][t] = s_error;
          stag[s][t] = &DataParser::parser_error;
        }
      after[s] = s_error;
      data [s] = &DataParser::white_spaces;
      etag [s] = &DataParser::end_tag;
    }


  // .....  <gnu-gama-data>  .........................................

  init(  s_start, t_gama_data, 
       //---------------------
       s_gama_data, 0, s_stop,
       &DataParser::gama_data, 0, 0);


  // .....  <g3-model>  ..............................................

  init(  s_gama_data, t_g3_model,
         //----------------------
         s_g3_model, 0, 0,
         &DataParser::g3_model, 0, &DataParser::g3_model);

  // .....  <g3-model> <vector>  ....................................

  init(  s_g3_model, t_vector,
         //-------------------
         s_g3_vector, s_g3_vector_after_czz, 0,
         &DataParser::g3_vector, 0, &DataParser::g3_vector);

  init(  s_g3_vector, t_from,
         //------------------
         s_g3_vector_from, 0, s_g3_vector_after_from,
         0, &DataParser::add_text, &DataParser::g3_vector_from);

  init(  s_g3_vector_after_from, t_to,
         //---------------------------
         s_g3_vector_to, 0, s_g3_vector_after_to,
         0, &DataParser::add_text, &DataParser::g3_vector_to);

  init(  s_g3_vector_after_from, t_to,
         //---------------------------
         s_g3_vector_to, 0, s_g3_vector_after_to,
         0, &DataParser::add_text, &DataParser::g3_vector_to);

  init(  s_g3_vector_after_to, t_dx,
         //-------------------------
         s_g3_vector_dx, 0, s_g3_vector_after_dx,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_dx, t_dy,
         //-------------------------
         s_g3_vector_dy, 0, s_g3_vector_after_dy,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_dy, t_dz,
         //-------------------------
         s_g3_vector_dz, 0, s_g3_vector_after_dz,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_dz, t_cxx,
         //--------------------------
         s_g3_vector_cxx, 0, s_g3_vector_after_cxx,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_cxx, t_cxy,
         //---------------------------
         s_g3_vector_cxy, 0, s_g3_vector_after_cxy,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_cxy, t_cxz,
         //---------------------------
         s_g3_vector_cxz, 0, s_g3_vector_after_cxz,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_cxz, t_cyy,
         //---------------------------
         s_g3_vector_cyy, 0, s_g3_vector_after_cyy,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_cyy, t_cyz,
         //---------------------------
         s_g3_vector_cyz, 0, s_g3_vector_after_cyz,
         0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_g3_vector_after_cyz, t_czz,
         //---------------------------
         s_g3_vector_czz, 0, s_g3_vector_after_czz,
         0, &DataParser::add_text, &DataParser::append_sp);

  // .....  <text>  ..................................................
 
  init(  s_gama_data, t_text,
       //--------------------
       s_text, 0, 0,
       0, &DataParser::add_text, &DataParser::text);
 
  // .....  <adj-input-data>  ........................................

  init(  s_gama_data, t_adj_input_data, 
       //------------------------------
       s_adj_input_data_1, s_adj_input_data_5, 0,
       &DataParser::adj_input_data, 0, &DataParser::adj_input_data,
       s_adj_input_data_4);

  // .....  <sparse-mat>  ............................................

  init(  s_adj_input_data_1, t_sparse_mat,
       //-------------------------------
       s_sparse_mat_1, s_sparse_mat_4, s_adj_input_data_2,
       0, 0, &DataParser::sparse_mat);

  init(  s_sparse_mat_1, t_rows, 
       //-----------------------
       s_sparse_mat_rows, 0, s_sparse_mat_2,
       0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_sparse_mat_2, t_cols,
       //-----------------------
       s_sparse_mat_cols, 0, s_sparse_mat_3,
       0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_sparse_mat_3, t_nonz,
       //-----------------------
       s_sparse_mat_nonz, 0, s_sparse_mat_4,
       0, &DataParser::add_text, &DataParser::sparse_mat_nonz);
 
  init(  s_sparse_mat_4, t_row, 
       //----------------------
       s_sparse_mat_row_1, s_sparse_mat_row_2, 0,
       &DataParser::sparse_mat_row, 0, &DataParser::sparse_mat_row);

  init(  s_sparse_mat_row_1, t_nonz,
       //--------------------------- 
       s_sparse_mat_row_nonz, 0, s_sparse_mat_row_2,
       0, &DataParser::add_text, &DataParser::sparse_mat_row_n);

  init(  s_sparse_mat_row_2, t_int, 
       //--------------------------
       s_sparse_mat_row_int, 0, s_sparse_mat_row_3,
       0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_sparse_mat_row_3, t_flt,
       //--------------------------
       s_sparse_mat_row_flt, 0, s_sparse_mat_row_2,
       0, &DataParser::add_text, &DataParser::sparse_mat_row_f);

  // ......  <block-diagonal>  .......................................

  init(  s_adj_input_data_2, t_block_diagonal, 
       //-----------------------------------
       s_block_diagonal_1, s_block_diagonal_3, s_adj_input_data_3,
       0, 0, &DataParser::block_diagonal);

  init(  s_block_diagonal_1, t_blocks,
       //-----------------------------
       s_block_diagonal_blocks, 0, s_block_diagonal_2,
       0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_block_diagonal_2, t_nonz,
       //---------------------------
       s_block_diagonal_nonz, 0, s_block_diagonal_3,
       0, &DataParser::add_text, &DataParser::block_diagonal_nonz);

  init(  s_block_diagonal_3, t_block,
       //----------------------------
       s_block_diagonal_block_1, s_block_diagonal_block_3, 0,
       0, 0, &DataParser::block_diagonal_block);

  init(  s_block_diagonal_block_1, t_dim,
       //--------------------------------
       s_block_diagonal_block_d, 0, s_block_diagonal_block_2,
       0, &DataParser::add_text, &DataParser::append_sp);

  init(  s_block_diagonal_block_2, t_width,
       //----------------------------------
       s_block_diagonal_block_w, 0, s_block_diagonal_block_3,
       0, &DataParser::add_text, &DataParser::block_diagonal_block_w);

  init(  s_block_diagonal_block_3, t_flt,  
       //--------------------------------
       s_block_diagonal_block_f, 0, 0,
       0, &DataParser::add_text, &DataParser::block_diagonal_vec_flt);

  // ......  <vector>  ...............................................

  init(  s_adj_input_data_3, t_vector, 
       //--------------------------
       s_vector_1, s_vector_2, s_adj_input_data_4,
       0, 0, &DataParser::vector);

  init(  s_vector_1, t_dim, 
       //------------------
       s_vector_dim, 0, s_vector_2,
       0, &DataParser::add_text, &DataParser::vector_dim);

  init(  s_vector_2, t_flt, 
       //------------------
       s_vector_flt, 0, 0,
       0, &DataParser::add_text, &DataParser::vector_flt);

  // ......  <array>  ................................................

  init(  s_adj_input_data_4, t_array,
       //--------------------------
       s_array_1, s_array_2, s_adj_input_data_5,
       0, 0, &DataParser::array);

  init(  s_array_1, t_dim,
       //-----------------
       s_array_dim, 0, s_array_2,
       0, &DataParser::add_text, &DataParser::array_dim);

  init(  s_array_2, t_int,
       //-----------------
       s_array_int, 0, 0,
       0, &DataParser::add_text, &DataParser::array_int);
         

  // .................................................................
}

// #######################################################
// #                                                     #
// # states:      s         n         z=n        a=s     #
// #              |         |         |          |       #
// # tag t:        <...t...>           </...t...>        #
// #                                                     #
// # functions:    [ Stag  ][  Data   ][  Etag  ]        #
// #                                                     #
// #######################################################

void DataParser::init(int s,   int t,           // current state, tag
                      int n,   int z,   int a,  // states: next , end, after
                      Stag s_, Data d_, Etag e_,
                      int z2)
{
  if (z == 0)  z = n;
  if (a == 0)  a = s;

  next [s][t] = n;
  after[z]    = a;
  
  if (s_) stag[s][t] = s_;
  else    stag[s][t] = &DataParser::start_tag;

  if (d_) data[n] = d_;

  if (e_) etag[z] = e_;

  if (z2)         // alternative end-state 
    { 
      after[z2] = a;
      etag [z2] = e_;
    }

}

DataParser::data_tag DataParser::tag(const char* c)
{
  switch (*c)
    {
    case 'a':
      if (!strcmp(c, "adj-input-data")) return t_adj_input_data;
      if (!strcmp(c, "array"         )) return t_array;
      break;
    case 'b':
      if (!strcmp(c, "block-diagonal")) return t_block_diagonal;
      if (!strcmp(c, "blocks"        )) return t_blocks;
      if (!strcmp(c, "block"         )) return t_block;
      break;
    case 'c' :
      if (!strcmp(c, "cols"          )) return t_cols;
      if (!strcmp(c, "cxx"           )) return t_cxx;
      if (!strcmp(c, "cxy"           )) return t_cxy;
      if (!strcmp(c, "cxz"           )) return t_cxz;
      if (!strcmp(c, "cyy"           )) return t_cyy;
      if (!strcmp(c, "cyz"           )) return t_cyz;
      if (!strcmp(c, "czz"           )) return t_czz;
      break;
    case 'd' :
      if (!strcmp(c, "dim"           )) return t_dim;
      if (!strcmp(c, "dx"            )) return t_dx;
      if (!strcmp(c, "dy"            )) return t_dy;
      if (!strcmp(c, "dz"            )) return t_dz;
      break;
    case 'f' :
      if (!strcmp(c, "flt"           )) return t_flt;
      if (!strcmp(c, "from"          )) return t_from;
      break;
    case 'g' :
      if (!strcmp(c, "g3-model"      )) return t_g3_model;
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
      if (!strcmp(c, "to"            )) return t_to;
      break;
    case 'v' :
      if (!strcmp(c, "vector"        )) return t_vector;
      break;
    case 'w' :
      if (!strcmp(c, "width"         )) return t_width;
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
      if (!isspace(s[len])) return error("### illegal text");
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
  objects.push_back( new DataObject::Text(text_buffer) );
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
  objects.push_back( new DataObject::AdjInput(data) );

  adj_sparse_mat = 0;
  adj_block_diagonal = 0;
  adj_array = 0;

  return end_tag(name);
}

// ......  <sparse-mat>  ...................................................

int DataParser::sparse_mat(const char *name)
{
  if (adj_sparse_mat && !adj_sparse_mat->check())
    error("### bad data in <sparse-mat>>");

  return end_tag(name);
}

int DataParser::sparse_mat_nonz(const char *name)
{
  std::size_t  rows, cols;
  istringstream inp(text_buffer.c_str());
  if (inp >> rows >> cols >> adj_sparse_mat_nonz)
    {
      text_buffer.erase();
      adj_sparse_mat = 
        new SparseMatrix<>(adj_sparse_mat_nonz, rows, cols);
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

// ......  <block-diagonal>  ...............................................

int DataParser::block_diagonal(const char *name)
{
  if (block_diagonal_blocks_)
    return error("### not enough <block> elements in <block-diagonal>");

  return end_tag(name);
}

int DataParser::block_diagonal_nonz(const char *name)
{
  istringstream inp(text_buffer.c_str());
  if (inp >> block_diagonal_blocks_ >> block_diagonal_nonz_)
    {
      text_buffer.erase();
      adj_block_diagonal = new BlockDiagonal<> 
        (block_diagonal_blocks_, block_diagonal_nonz_);
      return end_tag(name);
    }
  return error("### bad data in tags <blocks> / <nonz>");
}

int DataParser::block_diagonal_block_w(const char *name)
{
  istringstream inp(text_buffer.c_str());
  std::size_t dim, width;
  if ((inp >> dim >> width) && dim>0 && width>=0 && width<dim)
    {   
      block_diagonal_dim   = dim;
      block_diagonal_width = width;

      text_buffer.erase();
      bd_vector_dim = dim*(width+1) - width*(width+1)/2;
      bd_vector.reset(bd_vector_dim);
      bd_vector_iterator = bd_vector.begin();
      return end_tag(name);
    }
  return error("### bad data in tags <dim> / <width>");
}

int DataParser::block_diagonal_vec_flt(const char *name)
{
  if (bd_vector_dim == 0 || block_diagonal_nonz_ == 0)
    return error("### too many <flt> elements in <block-diagonal>");

  double flt;
  istringstream inp(text_buffer.c_str());
  if (inp >> flt)
    {
      bd_vector_dim--;
      block_diagonal_nonz_--;
      text_buffer.erase();
      *bd_vector_iterator++ = flt;
      
      return end_tag(name);
    }
  
  return error("### bad data format in a <flt> element in <block-diagonal>");
}

int DataParser::block_diagonal_block(const char *name)
{
  if (bd_vector_dim)
    return error("### not enough <flt> elements in <block-diagonal>");

  if (block_diagonal_blocks_ == 0)
    return error("### too many <block> elements in <block-diagonal>");

  block_diagonal_blocks_--;
  adj_block_diagonal->add_block(block_diagonal_dim, 
                                block_diagonal_width, bd_vector.begin());

  return end_tag(name);
}

// ......  <vector>  .......................................................

int DataParser::vector(const char *name)
{
  if (adj_vector_dim)
    return error("### not enough <flt> elements in <vector>");

  return end_tag(name);
}

int DataParser::vector_dim(const char *name)
{
  istringstream inp(text_buffer.c_str());
  if (inp >> adj_vector_dim)
    {
      text_buffer.erase();
      adj_vector.reset(adj_vector_dim);
      adj_vector_iterator = adj_vector.begin();
      
      return end_tag(name);
    }

  return error("### bad vector dimension in tag <dim>");
}

int DataParser::vector_flt(const char *name)
{
  if (adj_vector_dim == 0)
    return error("### too many <flt> elements in <vector>");

  double flt;
  istringstream inp(text_buffer.c_str());
  if (inp >> flt)
    {
      adj_vector_dim--;
      text_buffer.erase();
      *adj_vector_iterator++ = flt;
      
      return end_tag(name);
    }

  return error("### bad vector data in tag <flt>");
}

int DataParser::array(const char *name)
{
  if (adj_array_dim)
    return error("### not enough <int> elements in <array>");

  return end_tag(name);
}

int DataParser::array_dim(const char *name)
{
  istringstream inp(text_buffer.c_str());
  if (inp >> adj_array_dim)
    {
      text_buffer.erase();
      adj_array = new IntegerList<>(adj_array_dim);
      adj_array_iterator = adj_array->begin();
      return end_tag(name);
    }

  return error("### bad array dimension in tag <dim>");
}

int DataParser::array_int(const char *name)
{
  if (adj_array_dim == 0)
    return error("### too many <int> elements in <array>");

  int index;
  istringstream inp(text_buffer.c_str());
  if (inp >> index)
    {
      adj_array_dim--;
      text_buffer.erase();
      *adj_array_iterator++ = index;
      
      return end_tag(name);
    }

  return error("### bad array data in tag <int>");
}

int DataParser::g3_model(const char *name, const char **atts)
{
  no_attributes( name, atts );
  state = next[state][tag(name)];

  mg3 = new g3::Model;
  
  return 0;
}

int DataParser::g3_model(const char *name)
{
  objects.push_back( new DataObject::g3_model(mg3) );
  mg3 = 0;

  return  end_tag(name);
}

int DataParser::g3_vector(const char *name, const char **atts)
{
  no_attributes( name, atts );
  state = next[state][tag(name)];

  g3vec_from = "";
  g3vec_to   = "";

  return 0;
}

int DataParser::g3_vector(const char *name)
{
  stringstream istr(text_buffer);
  double dx, dy, dz, cxx, cxy, cxz, cyy, cyz, czz;

  if (!(istr >> dx >> dy >> dz >> cxx >> cxy >> cxz >> cyy >> cyz >> czz))
    {
      return error("### bad format of numerical data in <vector>");
    }

  using namespace g3;
  
  cout << dx<< " " <<  dy<< " " <<  dz<< " " 
       <<  cxx<< " " <<  cxy<< " " <<  cxz<< " " 
       <<  cyy<< " " <<  cyz<< " " <<  czz << endl;

  Vector* v = new Vector(dx, dy, dz);
  v->name[0] = g3vec_from;
  v->name[1] = g3vec_to;

  Model::ObservationData *obs = mg3->obs;
  Vectors* vectors = new Vectors(obs);

  vectors->add(v);
  obs->CL.push_back(vectors);

  return  end_tag(name);
}

int DataParser::g3_vector_from(const char *name)
{
  g3vec_from = text_buffer;
  text_buffer.erase();

  return  end_tag(name);
}

int DataParser::g3_vector_to(const char *name)
{
  g3vec_to = text_buffer;
  text_buffer.erase();

  return  end_tag(name);
}

#endif
