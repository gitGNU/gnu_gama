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
 *  $Id: dataparser.h,v 1.6 2003/05/15 18:53:15 cepek Exp $
 */

#ifndef GNU_Gama_GaMa_XML_DataParser__data_parser__dataparser___h_
#define GNU_Gama_GaMa_XML_DataParser__data_parser__dataparser___h_

#include <gnu_gama/xml/baseparser.h>
#include <gnu_gama/xml/dataobject.h>
#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/list.h>
#include <gnu_gama/exception.h>
#include <cstddef>
#include <string>
#include <list>

namespace GNU_gama {
  
  class DataParser : public BaseParser<Exception::parser>
    {
    public:
      
      DataParser(List<DataObject::Base*>&);
      ~DataParser();
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
      
      List<DataObject::Base*>& objects;
      
      enum parser_state 
        {
          s_error,
          s_start,
          s_gama_data,
          s_g3_model,
          s_text,
          s_adj_input_data_1,
          s_adj_input_data_2,
          s_adj_input_data_3,
          s_adj_input_data_4,
          s_adj_input_data_5,
          s_sparse_mat_1,
          s_sparse_mat_rows,
          s_sparse_mat_2,
          s_sparse_mat_cols,
          s_sparse_mat_3,
          s_sparse_mat_nonz,
          s_sparse_mat_4,
          s_sparse_mat_row_1,
          s_sparse_mat_row_nonz,
          s_sparse_mat_row_2,
          s_sparse_mat_row_int,
          s_sparse_mat_row_3,
          s_sparse_mat_row_flt,
          s_block_diagonal_1,
          s_block_diagonal_blocks,
          s_block_diagonal_2,
          s_block_diagonal_nonz,
          s_block_diagonal_3,
          s_block_diagonal_block_1,
          s_block_diagonal_block_d,
          s_block_diagonal_block_2,
          s_block_diagonal_block_w,
          s_block_diagonal_block_3,
          s_block_diagonal_block_f,
          s_vector_1,
          s_vector_dim,
          s_vector_2,
          s_vector_flt,
          s_array_1,
          s_array_dim,
          s_array_2,
          s_array_int,
          s_stop       
        }; 
      
      enum data_tag 
        {
          t_adj_input_data,
          t_array,
          t_block,
          t_block_diagonal,
          t_blocks,
          t_cols,
          t_dim,
          t_flt,
          t_g3_model,
          t_gama_data,
          t_int,
          t_nonz,
          t_rows,
          t_row,
          t_sparse_mat,
          t_text,
          t_vector,
          t_width,
          t_unknown
        };
      
      data_tag tag(const char *name);
      
      typedef int (DataParser::*Stag)(const char *name, const char **atts);
      typedef int (DataParser::*Data)(const char *name, int len);
      typedef int (DataParser::*Etag)(const char *name);

      Stag stag[s_stop+1][t_unknown+1];
      Data data[s_stop+1];                 
      Etag etag[s_stop+1];

      int next [s_stop+1][t_unknown+1];
      int after[s_stop+1]; 

      int gama_data             (const char *name, const char **atts);
      int g3_model              (const char *name, const char **atts);
      int g3_model              (const char *name);
      int text                  (const char *name);
      int adj_input_data        (const char *name, const char **atts);
      int adj_input_data        (const char *name);
      int sparse_mat            (const char *name);
      int sparse_mat_nonz       (const char *name);
      int sparse_mat_row        (const char *name, const char **atts);
      int sparse_mat_row        (const char *name);
      int sparse_mat_row_n      (const char *name);
      int sparse_mat_row_f      (const char *name);
      int block_diagonal        (const char *name);
      int block_diagonal_nonz   (const char *name);
      int block_diagonal_block_w(const char *name);
      int block_diagonal_vec_flt(const char *name);
      int block_diagonal_block  (const char *name);
      int vector                (const char *name);
      int vector_dim            (const char *name);
      int vector_flt            (const char *name);
      int array                 (const char *name);
      int array_dim             (const char *name);
      int array_int             (const char *name);
      
      int add_text     (const char *name, int len);
      int end_tag      (const char *name);
      int no_attributes(const char *name, const char **atts);
      int parser_error (const char *name, const char **atts);
      int start_tag    (const char *name, const char **atts);
      int white_spaces (const char *name, int len);
      int append_sp    (const char *name);

      void init(int state, int tag, 
                int next_state, int end_state, int after_state,
                Stag, Data, Etag,
                int end_state2=0);


      // DataObject::g3_model

      g3::Model*       mg3;

      // DataObject::Text

      std::string      text_buffer;

      // DataObject::AdjInput

      GNU_gama::SparseMatrix <> *adj_sparse_mat;
      GNU_gama::BlockDiagonal<> *adj_block_diagonal;
      Vec              adj_vector;  
      Vec::iterator    adj_vector_iterator;  
      std::size_t      adj_vector_dim;
      GNU_gama::IntegerList<>   *adj_array;
      GNU_gama::IntegerList<>::iterator adj_array_iterator;
      std::size_t      adj_array_dim;
      std::size_t      adj_sparse_mat_nonz;
      std::size_t      adj_sparse_mat_row_nonz;
      std::size_t      block_diagonal_blocks_;
      std::size_t      block_diagonal_nonz_;  
      std::size_t      block_diagonal_dim;
      std::size_t      block_diagonal_width;
      Vec              bd_vector;
      Vec::iterator    bd_vector_iterator;  
      std::size_t      bd_vector_dim;

    };
}

#endif
