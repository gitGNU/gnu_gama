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
 *  $Id: dataparser.h,v 1.2 2002/10/18 20:52:29 cepek Exp $
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
    ~DataParser();
     
    int characterDataHandler(const char* s, int len);
    int startElement(const char *cname, const char **atts);
    int endElement(const char * name);
    
  private: 

    std::list<DataObject*>& objects;

    enum parser_state 
      {
        state_error,
        state_start,
        state_gama_data,
        state_stop       
      } next[state_stop]; 

    enum data_tag 
      {
        tag_gama_data,
        tag_text,
        tag_unknown
      };

    typedef int (DataParser::*FUN)(const char *cname, const char **atts);
    FUN fun[state_stop+1][tag_unknown+1];

    data_tag tag(const char* cname);

    int t_error    (const char *cname, const char **atts);
    int t_gama_data(const char *cname, const char **atts);

  };  // class DataParser
}     // namespace GaMaLib


#endif

















