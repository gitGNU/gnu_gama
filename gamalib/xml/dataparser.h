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
 *  $Id: dataparser.h,v 1.8 2003/01/06 17:44:12 cepek Exp $
 */

#ifndef GaMaLib_GaMa_XML_DataParser__data_parser__dataparser___h_
#define GaMaLib_GaMa_XML_DataParser__data_parser__dataparser___h_

#include <gamalib/xml/baseparser.h>
#include <gamalib/xml/dataobject.h>
#include <gamalib/local/gamadata.h>
#include <cstddef>
#include <string>
#include <list>

namespace GaMaLib {
  
  class DataParser : public BaseParser
    {
    public:
      
      DataParser(std::list<DataObject*>&);
      ~DataParser()
        {
        } 
      int startElement(const char *name, const char **atts)
        {
          return (this->*stag[state][tag(name)])(name, atts);
        }  
      int characterDataHandler(const char *s, int len)
        { 
          return (this->*data[state])(s, len);
        } 
      int endElement(const char *name)
        {
          return (this->*etag[state])(name);
        }
      
    private: 
      
      std::list<DataObject*>& objects;
      
      enum parser_state 
        {
          s_error,
          s_start,
          s_gama_data,
          s_text,
          s_adj_input_data,
          s_sparse_mat_1,
          s_sparse_mat_rows,
          s_sparse_mat_2,
          s_sparse_mat_cols,
          s_sparse_mat_3,
          s_sparse_mat_nonz,
          s_sparse_mat_4,
          s_sparse_mat_row_1,
          s_sparse_mat_row_n,
          s_sparse_mat_row_2,
          s_sparse_mat_row_i,
          s_sparse_mat_row_3,
          s_sparse_mat_row_f,
          s_stop       
        }; 
      
      enum data_tag 
        {
          t_adj_input_data,
          t_cols,
          t_flt,
          t_gama_data,
          t_int,
          t_nonz,
          t_rows,
          t_row,
          t_sparse_mat,
          t_text,
          t_unknown
        };
      
      data_tag tag(const char* name);
      
      typedef int (DataParser::*Stag)(const char *name, const char **atts);
      typedef int (DataParser::*Data)(const char *s, int len);
      typedef int (DataParser::*Etag)(const char *name);

      Stag stag[s_stop+1][t_unknown+1];
      Data data[s_stop+1];                 
      Etag etag[s_stop+1];

      int next [s_stop+1][t_unknown+1];
      int after[s_stop+1]; 

      int gama_data       (const char *name, const char **atts);
      int text            (const char* s);
      int adj_input_data  (const char *name, const char **atts);
      int adj_input_data  (const char *name);
      int sparse_mat      (const char *name);
      int sparse_mat_nonz (const char *name);
      int sparse_mat_row  (const char *name, const char **atts);
      int sparse_mat_row  (const char *name);
      int sparse_mat_row_n(const char *name);
      int sparse_mat_row_f(const char *name);
      

      int add_text     (const char* s, int len);
      int end_tag      (const char *name);
      int no_attributes(const char *name, const char **atts);
      int parser_error (const char *name, const char **atts);
      int start_tag    (const char *name, const char **atts);
      int white_spaces (const char* s, int len);
      int append_sp    (const char *name);

      
      
      std::string      text_buffer;
      SparseMatrix <> *adj_sparse_mat;
      BlockDiagonal<> *adj_block_diagonal;
      Vec              adj_vector;  
      IntegerList<>   *adj_array;
      std::size_t      adj_sparse_mat_nonz;
      std::size_t      adj_sparse_mat_row_nonz;
      
    };
}

#endif
