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
 *  $Id: baseparser.h,v 1.1 2003/05/10 13:13:36 cepek Exp $
 */

#ifndef GaMaLib_GaMa__XML__BASE_Base_base__PARSER_Parser_parser__h_
#define GaMaLib_GaMa__XML__BASE_Base_base__PARSER_Parser_parser__h_


// BaseParser is just a simple C++ wrapper for XML parser expat

#include <expat/xmlparse/xmlparse.h>

#include <gamalib/xml/dataobject.h>
#include <gamalib/local/gamadata.h>
#include <string>
#include <list>

namespace GNU_gama {

  class ParserException : public GaMaLib::Exception 
  {
  public:

    int line, error_code;

    ParserException(std::string s, int r, int c)
      : GaMaLib::Exception(s), line(r), error_code(c) 
      {
      }

  };
  
  class BaseParser 
  {
  public:
    
    BaseParser();
    virtual ~BaseParser();
    
    // expat parser interface
    
    void xml_parse(const char *s, int len, int  isFinal); 
    void xml_parse(const std::string& s, bool isFinal) 
    {
      xml_parse(s.c_str(), s.length(), isFinal ? 1 : 0);
    }
    virtual int characterDataHandler(const char* s, int len) = 0;
    virtual int startElement(const char *cname, const char **atts) = 0;
    virtual int endElement(const char * name) = 0;

  protected: 
    
    XML_Parser  parser;
    int         state;      /*  state_error must be 0  */
    
    int error(const char* text);
    int error(const std::string& s)   { return error(s.c_str()); }

    bool toDouble(const std::string&, double&     ) const;
    bool toIndex (const std::string&, std::size_t&) const;
      
  private:

    std::string errString;
    int         errLineNumber;  
    int         errCode;              // -1 bad data in gkf; 0 OK; >0 expat

  };  // class DataParser
}     // namespace GaMaLib


#endif

















