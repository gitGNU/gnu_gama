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
 *  $Id: g3_obs_dist.cpp,v 1.1 2003/03/28 22:07:31 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation.h>
#include <gnu_gama/g3/g3_model.h>
#include <cmath>

using namespace GNU_gama::g3;
using namespace std;


double Distance::parlist_value() const
{
  Parameter** p = parlist.begin();
  double b1 = (*p++)->value(time);
  double l1 = (*p++)->value(time);
  double h1 = (*p++)->value(time);
  double b2 = (*p++)->value(time);
  double l2 = (*p++)->value(time);
  double h2 = (*p++)->value(time);

  double x1, y1, z1, x2, y2, z2;
  ellipsoid->blh2xyz(b1, l1, h1, x1, y1, z1);
  ellipsoid->blh2xyz(b2, l2, h2, x2, y2, z2);
  x2 -= x1;
  y2 -= y2;
  z2 -= z1;

  return sqrt(x2*x2 + y2*y2 + z2*z2);
}


void Distance::parlist_init(Model* model)
{
  ellipsoid = &model->ellipsoid;

  if (!active())  return;
 
  Point* from = model->points.find(name[0]);
  Point* to   = model->points.find(name[1]);

  if (from == 0 || to == 0)
    {
      set_active(false);
      return;
    }


  if (from->B->active() && from->L->active() && from->H->active() && 
      to->  B->active() && to->  L->active() && to->  H->active()  )
    {
      Parameter** b = parlist.begin();
      *b++ = from->B;
      *b++ = from->L;
      *b++ = from->H;
      *b++ = to->B;
      *b++ = to->L;
      *b++ = to->H;
    }
  else
    {
      set_active(false);
      return;
    }
}


double Distance::derivative(Parameter* p)
{
  DistanceAnalyticalDerivative* ad = 
    dynamic_cast<DistanceAnalyticalDerivative*>(p);

  if (ad)
    return ad->analytical_derivative(this);
  else
    return numerical_derivative(p);
}
