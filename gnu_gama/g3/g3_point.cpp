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
 *  $Id: g3_point.cpp,v 1.15 2003/04/11 09:38:26 cepek Exp $
 */

#include <gnu_gama/g3/g3_point.h>
#include <gnu_gama/g3/g3_observation.h>
#include <gnu_gama/g3/g3_model.h>
#include <cmath>

using namespace GNU_gama::g3;


Point::~Point()
{
  delete N;
  delete E;
  delete U;
}


Point::Point() : parlist(3)
{
  N = new Parameter_N(this);
  E = new Parameter_E(this);
  U = new Parameter_U(this);

  set_unused();

  Parameter** p = parlist.begin();
  *p++ = N;
  *p++ = E;
  *p++ = U;
}


Point::Point(const Point& point) : parlist(3)
{
  name   = point.name;
  common = point.common;

  N = new Parameter_N(*point.N);
  E = new Parameter_E(*point.E);
  U = new Parameter_U(*point.U);

  Parameter** p = parlist.begin();
  *p++ = N;
  *p++ = E;
  *p++ = U;
}


Point& Point::operator=(const Point& point)
{
  if (this != &point)
    {
      delete N;
      delete E;
      delete U;

      name   = point.name;
      common = point.common;
      
      N = new Parameter_N(*point.N);
      E = new Parameter_E(*point.E);
      U = new Parameter_U(*point.U);
      
      N->set_point(this);
      E->set_point(this);
      U->set_point(this);

      Parameter** p = parlist.begin();
      *p++ = N;
      *p++ = E;
      *p++ = U;
    }

  return *this;
}


void Point::set_unused()
{
  N->set_unused();
  E->set_unused();
  U->set_unused();
}

void Point::set_fixed_horizontal_position()
{
  N->set_fixed();
  E->set_fixed();
}

void Point::set_fixed_height()
{
  U->set_fixed();
}

void Point::set_fixed_position()
{
  N->set_fixed();
  E->set_fixed();
  U->set_fixed();
}

void Point::set_free_horizontal_position()
{
  N->set_free();
  E->set_free();
}

void Point::set_free_height()
{
  U->set_free();
}

void Point::set_free_position()
{
  N->set_free();
  E->set_free();
  U->set_free();
}

void Point::set_constr_horizontal_position()
{
  N->set_constr();
  E->set_constr();
}

void Point::set_constr_height()
{
  U->set_constr();
}

void Point::set_constr_position()
{
  N->set_constr();
  E->set_constr();
  U->set_constr();
}

bool Point::unused() const
{
  return N->unused() && E->unused() && U->unused();
}

bool Point::fixed_horizontal_position() const
{
  return N->fixed() && E->fixed();
}

bool Point::fixed_height() const
{
  return U->fixed();
}

bool Point::fixed_position() const
{
  return N->fixed() && E->fixed() && U->fixed();
}

bool Point::free_horizontal_position() const
{
  return N->free() && E->free();
}

bool Point::free_height() const
{
  return U->free();
}

bool Point::free_position() const
{
  return N->free() && E->free() && U->free();
}

bool Point::constr_horizontal_position() const
{
  return N->constr() && E->constr();
}

bool Point::constr_height() const
{
  return U->constr();
}

bool Point::constr_position() const
{
  return N->constr() && E->constr() && U->constr();
}

void Point::set_blh(double b, double l, double h)
{
  B.set_init_value(b);
  L.set_init_value(l);
  H.set_init_value(h);
  
  double x, y, z;
  common->ellipsoid.blh2xyz(b, l, h, x, y, z);
  X.set_init_value(x);
  Y.set_init_value(y);
  Z.set_init_value(z);

  transformation_matrix(b, l);
}

void Point::set_xyz(double x, double y, double z)
{
  X.set_init_value(x);
  Y.set_init_value(y);
  Z.set_init_value(z);
  
  double b, l, h;;
  common->ellipsoid.xyz2blh(x, y, z, b, l, h);
  B.set_init_value(b);
  L.set_init_value(l);
  H.set_init_value(h);

  transformation_matrix(b, l);
}

void Point::transformation_matrix(double b, double l)
{
  using namespace std;

  r11 = -sin(l);
  r12 = -sin(b)*cos(l);
  r13 =  cos(b)*cos(l);
  
  r21 =  cos(l);
  r22 = -sin(b)*sin(l);
  r23 =  cos(b)*sin(l);

  r31 = 0.0;
  r32 = cos(b);
  r33 = sin(b);
}


double Point::x_transform(double n, double e, double u)
{
  return r11*n + r12*e + r13*u;
}

double Point::y_transform(double n, double e, double u)
{
  return r21*n + r22*e + r23*u;
}

double Point::z_transform(double n, double e, double u)
{
  return r31*n + r32*e + r33*u;
}


void Point::set_diff_XYZ(double dx, double dy, double dz)
{
  dX = dx;  dY = dy;  dZ = dz;
}

double Point::diff_N() const
{
  return r11*dX + r21*dY + r31*dZ;
}

double Point::diff_E() const
{
  return r12*dX + r22*dY + r32*dZ;
}

double Point::diff_U() const
{
  return r13*dX + r23*dY + r33*dZ;
}


// ----------------------------------------------------------------------

double Parameter_N::derivative(Distance*)
{
  return point()->diff_N();
}

double Parameter_E::derivative(Distance*)
{
  return point()->diff_E();
}

double Parameter_U::derivative(Distance*)
{
  return point()->diff_U();
}

double Parameter_U::derivative(HeightDiff* hd)
{
  if (!free()) return 0;

  Parameter**  p    = hd->parlist.begin();  
  Parameter_U* from = static_cast<Parameter_U*>(*p);

  return this == from ? -1 : +1;
}

