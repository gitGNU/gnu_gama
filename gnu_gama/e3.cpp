/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2004  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: e3.cpp,v 1.1 2004/05/13 17:50:27 cepek Exp $
 */

#include <gnu_gama/e3.h>
#include <cmath>

using namespace GNU_gama;

void E3::operator += (const E3& v)
{
  e1 += v.e1;  e2 += v.e2; e3 += v.e3;
}
void E3::operator -= (const E3& v)
{
  e1 -= v.e1;  e2 -= v.e2; e3 -= v.e3;
}
void E3::operator *= (double d)
{
  e1 *= d;  e2 *= d;  e3 *= d;
}
void E3::set(double a, double b, double c)
{
  e1 = a;  e2 = b;  e3 = c;
}
void E3::add(double a, double b, double c)
{
  e1 += a;  e2 += b;  e3 += c;
}
void E3::sub(double a, double b, double c)
{
  e1 -= a;  e2 -= b;  e3 -= c;
}
double E3::dot(const E3& v) const
{
  return e1*v.e1 + e2*v.e2 + e3*v.e3;
}
void E3::cross(const E3& a, const E3& b)
{
  set(+a.e2*b.e3 - b.e2*a.e3,
      -a.e1*b.e3 + b.e1*a.e3,
      +a.e1*b.e2 - b.e1*a.e2);
}

double GNU_gama::angle(const GNU_gama::E3& a, const GNU_gama::E3& b)
{
  E3 x;  x.cross(a, b);

  const double s = x.dot(x);
  const double c = a.dot(b);

  return std::atan2(s, c);
}

