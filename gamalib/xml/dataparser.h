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
 *  $Id: dataparser.h,v 1.7 2003/01/05 18:02:33 cepek Exp $
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
      int characterDataHandler(const char* s, int len)
        {
          return (this->*data[state])(s, len);
        }
      int startElement(const char *name, const char **atts)
        {
          data_tag t = tag(name);
          int s = state;
          stack.push(state);
          state = next[state][t];
          return (this->*func[s][t])(name, atts);
        }
      int endElement(const char * /*name*/)
        {
          if (Ende e = ende[state]) (this->*e)();
          if (state) state = stack.pop();
          return state;
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
          s_sparse_mat,
          s_sparse_mat_rows,
          s_sparse_mat_cols,
          s_sparse_mat_nonz,
          s_sparse_mat_row,
          s_sparse_mat_row_n,
          s_sparse_mat_row_i,
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
      
      data_tag tag(const char* cname);
      
      class Stack {
        int N;
        int buf[s_stop];
      public:
        Stack() : N(0) {}
        void push(int s) { buf[N++] = s;    }
        int  pop ()      { return buf[--N]; }
      };
      Stack stack;
      
      typedef int (DataParser::*Func)(const char *cname, const char **atts);
      Func func[s_stop+1][t_unknown+1];
      int  next[s_stop+1][t_unknown+1];
      
      typedef int (DataParser::*Data)(const char* s, int len);
      Data data[s_stop+1];
      
      typedef int (DataParser::*Ende)();
      Ende ende[s_stop+1];
      
      int parser_error (const char *cname, const char **atts);
      int no_attributes(const char *cname, const char **atts);
      int white_spaces (const char* s, int len);
      int add_text     (const char* s, int len);
      int add_space    ();
      
      int gama_data(const char *cname, const char **atts);
      
      int text(const char *cname, const char **atts);
      int text(const char* s, int len);
      int text();
      
      int adj_input_data(const char *cname, const char **atts);
      int adj_input_data();
      
      int sparse_mat      ();
      int sparse_mat_nonz ();
      int sparse_mat_row  (const char *cname, const char **atts);
      int sparse_mat_row_n();
      int sparse_mat_row_f();
      
      std::string      text_buffer;
      SparseMatrix <> *tmp_sparse_mat;
      BlockDiagonal<> *tmp_block_diagonal;
      Vec              tmp_vector;  
      IntegerList<>   *tmp_array;
      std::size_t      tmp_sparse_mat_nonz;
      std::size_t      tmp_sparse_mat_row_nonz;
      
    };
}

#endif
