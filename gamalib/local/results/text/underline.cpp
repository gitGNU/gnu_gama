/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2001  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: underline.cpp,v 1.4 2003/05/10 13:43:03 cepek Exp $
 */

#include <string>
#include <gamalib/local/results/text/underline.h>
#include <gnu_gama/xml/encoding.h>


namespace GaMaLib {

std::string underline(std::string text, char c)
{
  int i;
  std::string s;
  unsigned char* p = (unsigned char*)text.c_str();
  while (*p)
    {
      p += GNU_gama::Utf8Decode(i, p);
      s += c;
    }
  return s;
}

std::string set_width(std::string s, int n)
{
  int N=0, i;
  unsigned char* p = (unsigned char*)s.c_str();
  while (*p)
    {
      p += GNU_gama::Utf8Decode(i, p);
      N++;
    }
  std::string t(s);
  while (n-- > N) t += ' ';
  return t;
}

}
