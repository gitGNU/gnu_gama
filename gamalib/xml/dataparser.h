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
 *  $Id: dataparser.h,v 1.5 2003/01/04 22:11:48 cepek Exp $
 */

#ifndef GaMaLib_GaMa_XML_Data_Object___parser__h_
#define GaMaLib_GaMa_XML_Data_Object___parser__h_

#include <gamalib/xml/baseparser.h>
#include <gamalib/xml/dataobject.h>
#include <gamalib/local/gamadata.h>
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
        return (this->*dat[state])(s, len);
      }
    int startElement(const char *name, const char **atts)
      {
        data_tag t = tag(name);
        int s = state;
        stack.push(state);
        state = next[state][t];
        return (this->*fun[s][t])(name, atts);
      }
    int endElement(const char * /*name*/)
      {
        if (state) state = stack.pop();
        return state;
      }

  private: 

    std::list<DataObject*>& objects;

    enum parser_state 
      {
        state_error,
        state_start,
        state_gama_data,
        state_text,
        state_adj_input_data,
        state_adj_input_data_sm1,
        state_adj_input_data_sm1_rows,
        state_adj_input_data_sm2,
        state_adj_input_data_sm2_cols,
        state_adj_input_data_sm3,
        state_adj_input_data_sm3_nonz,
        state_adj_input_data_sm,
        state_stop       
      }; 

    enum data_tag 
      {
        tag_adj_input_data,
        tag_gama_data,
        tag_sparse_mat,
        tag_text,
        tag_unknown
      };

    data_tag tag(const char* cname);

    class Stack {
      int N;
      int buf[state_stop];
    public:
      Stack() : N(0) {}
      void push(int s) { buf[N++] = s;    }
      int  pop ()      { return buf[--N]; }
    };
    Stack stack;
    
    typedef int (DataParser::*FUN)(const char *cname, const char **atts);
    FUN fun [state_stop+1][tag_unknown+1];
    int next[state_stop+1][tag_unknown+1];

    int t_error         (const char *cname, const char **atts);
    int t_no_attributes (const char *cname, const char **atts);
    int t_gama_data     (const char *cname, const char **atts);
    int t_text          (const char *cname, const char **atts);
    int t_adj_input_data(const char *cname, const char **atts);

    AdjInputDataObject *ptr_adj_input_data;

    typedef int (DataParser::*DATA)(const char* s, int len);
    DATA dat[state_stop+1];

    int d_ws  (const char* s, int len);
    int d_text(const char* s, int len);
  };
}

#endif

















