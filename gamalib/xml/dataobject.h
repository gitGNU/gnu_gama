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
 *  $Id: dataobject.h,v 1.5 2003/01/04 15:51:51 cepek Exp $
 */

#ifndef GaMaLib_GaMa_XML_Data_Object__object___h_
#define GaMaLib_GaMa_XML_Data_Object__object___h_

#include <string>
#include <sstream>
#include <gamalib/adj/adj.h>

namespace GaMaLib {

  class DataObject {
  public:

    virtual ~DataObject() 
      {
      }
    virtual std::string xml() const = 0;

    static  std::string xml_begin();
    static  std::string xml_end();
  };


  class TextDataObject : public DataObject {
  public:
  
    std::string text;
  
    TextDataObject() 
      {
      }    
    TextDataObject(std::string s) : text(s) 
      {
      }    
    std::string xml() const 
      {
        if (!text.empty())
          return "\n<text>" + text + "</text>\n";

        return "";
      }
  };


  class AdjInputDataObject : public DataObject {
  public:
  
    AdjInputData *data;
  
    AdjInputDataObject() : data(0)
      {
      }    
    AdjInputDataObject(AdjInputData *d) : data(d)
      {
      }    
    std::string xml() const 
      {
        if (data) 
          {
            std::stringstream out;
            data->write_xml(out);
            
            return out.str();
          }

        return "";
      }
  };

}       // namespace GaMaLib


#endif

















