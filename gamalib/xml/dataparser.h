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
 *  $Id: dataparser.h,v 1.1 2002/10/17 17:24:55 cepek Exp $
 */

#ifndef GaMaLib_GaMa_XML_Data_Object___parser__h_
#define GaMaLib_GaMa_XML_Data_Object___parser__h_


// DataParser is just a simple C++ wrapper for XML parser expat

#include <gamalib/xml/baseparser.h>
#include <gamalib/xml/dataobject.h>
#include <gamalib/local/gamadata.h>
#include <string>
#include <list>

namespace GaMaLib {

  class DataParserException : public GaMaLib::Exception 
  {
  public:
    int line;
    DataParserException(std::string s, int r) 
      : GaMaLib::Exception(s), line(r) {}
  };
  
  class DataParser 
  {
  public:
    
    std::string errString;
    int         errLineNumber;  
    int         errCode;              // -1 bad data in gkf; 0 OK; >0 expat
    
    // constructor and destructor
    
    DataParser(std::list<DataObject*>&);
    ~DataParser();
    
    // expat parser interface
    
    void xml_parse(const std::string& s, bool isFinal) 
    {
      xml_parse(s.c_str(), s.length(), isFinal ? 1 : 0);
    }
    void xml_parse(const char *s, int len, int  isFinal) 
    { 
      int err = XML_Parse(parser, s, len, isFinal);
      if (err == 0)
        {
          // fatal error
               
          errString=std::string(XML_ErrorString(XML_GetErrorCode(parser)));
          errCode  =XML_GetErrorCode(parser);
          errLineNumber = XML_GetCurrentLineNumber(parser);
          
          throw DataParserException(errString, errLineNumber);
        }
      
      if (state == state_error)
        {
          errCode = -1;   
          // errLineNumber is set by function  error("...");    
          throw DataParserException(errString, errLineNumber);
        }
    }
    
    // int gkf_characterDataHandler(const char* s, int len);
    // int gkf_startElement(const char *cname, const char **atts);
    // int gkf_endElement(const char * name);
    
  private: 

    XML_Parser  parser;

    enum parser_state {
      state_error,
      state_start,
      state_stop
    } state;

    int error(std::string s) { return error(s.c_str()); }
    int error(const char* text)
    {
      // store only the first detected error
      if(errCode) return 1;
      
      errString = std::string(text);
      errCode   = -1;
      errLineNumber = XML_GetCurrentLineNumber(parser);
      state = state_error;
      return 1;
    }
    
    bool toDouble(const std::string&, double&) const;
    bool toIndex (const std::string&, Index& ) const;
      

  };  // class DataParser
}     // namespace GaMaLib


#endif

















