/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
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
 *  $Id: g3_point.cpp,v 1.10 2003/03/28 22:07:31 cepek Exp $
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
  B = new Parameter_B;
  L = new Parameter_L;
  H = new Parameter_H;

  set_unused();

  Parameter** p = parlist.begin();
  *p++ = B;
  *p++ = L;
  *p++ = H;
}


Point::Point(const Point& point) : parlist(3)
{
  name   = point.name;
  common = point.common;

  B = new Parameter_B(*point.B);
  L = new Parameter_L(*point.L);
  H = new Parameter_H(*point.H);

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
      
      B = new Parameter_B(*point.B);
      L = new Parameter_L(*point.L);
      H = new Parameter_H(*point.H);
      
      Parameter** p = parlist.begin();
      *p++ = B;
      *p++ = L;
      *p++ = H;
    }

  return *this;
}


void Point::set_unused()
{
  B->set_unused();
  L->set_unused();
  H->set_unused();
}

void Point::set_fixed_horizontal_position()
{
  B->set_fixed();
  L->set_fixed();
}

void Point::set_fixed_height()
{
  H->set_fixed();
}

void Point::set_fixed_position()
{
  B->set_fixed();
  L->set_fixed();
  H->set_fixed();
}

void Point::set_free_horizontal_position()
{
  B->set_free();
  L->set_free();
}

void Point::set_free_height()
{
  H->set_free();
}

void Point::set_free_position()
{
  B->set_free();
  L->set_free();
  H->set_free();
}

void Point::set_constr_horizontal_position()
{
  B->set_constr();
  L->set_constr();
}

void Point::set_constr_height()
{
  H->set_constr();
}

void Point::set_constr_position()
{
  B->set_constr();
  L->set_constr();
  H->set_constr();
}

bool Point::unused() const
{
  return B->unused() && L->unused() && H->unused();
}

bool Point::fixed_horizontal_position() const
{
  return B->fixed() && L->fixed();
}

bool Point::fixed_height() const
{
  return H->fixed();
}

bool Point::fixed_position() const
{
  return B->fixed() && L->fixed() && H->fixed();
}

bool Point::free_horizontal_position() const
{
  return B->free() && L->free();
}

bool Point::free_height() const
{
  return H->free();
}

bool Point::free_position() const
{
  return B->free() && L->free() && H->free();
}

bool Point::constr_horizontal_position() const
{
  return B->constr() && L->constr();
}

bool Point::constr_height() const
{
  return H->constr();
}

bool Point::constr_position() const
{
  return B->constr() && L->constr() && H->constr();
}


