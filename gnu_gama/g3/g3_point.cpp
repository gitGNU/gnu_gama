/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama library.
    
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
 *  $Id: g3_point.cpp,v 1.6 2003/03/22 21:55:29 cepek Exp $
 */

#include <gnu_gama/g3/g3_point.h>


using namespace GNU_gama::g3;


Point::~Point()
{
  delete B;
  delete L;
  delete H;
}


Point::Point() : parlist(3)
{
  B = new Parameter_position;
  L = new Parameter_position;
  H = new Parameter_height;

  Parameter** p = parlist.begin();
  *p++ = B;
  *p++ = L;
  *p++ = H;
}


Point::Point(const Point& point) : parlist(3)
{
  name   = point.name;
  common = point.common;

  B = new Parameter_position;
  L = new Parameter_position;
  H = new Parameter_height;

  Parameter** p = parlist.begin();
  *p++ = B;
  *p++ = L;
  *p++ = H;
}




Point& Point::operator=(const Point& point)
{
  if (this != &point)
    {
      delete B;
      delete L;
      delete H;

      name   = point.name;
      common = point.common;
      
      B = new Parameter_position;
      L = new Parameter_position;
      H = new Parameter_height;
      
      Parameter** p = parlist.begin();
      *p++ = B;
      *p++ = L;
      *p++ = H;
    }

  return *this;
}


