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
 *  $Id: g3_obs_dist.cpp,v 1.2 2003/04/08 16:41:51 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation.h>
#include <gnu_gama/g3/g3_model.h>
#include <cmath>

using namespace GNU_gama::g3;
using namespace std;


double Distance::parlist_value() const
{
  Parameter** p = parlist.begin();

  Parameter_N* N_from = static_cast<Parameter_N*>(*p);

  double n1 = (*p++)->value(time);
  double e1 = (*p++)->value(time);
  double u1 = (*p++)->value(time);

  Parameter_N* N_to   = static_cast<Parameter_N*>(*p);

  double n2 = (*p++)->value(time);
  double e2 = (*p++)->value(time);
  double u2 = (*p++)->value(time);

  Point* from = N_from->point();
  Point* to   = N_to  ->point();

  double x1 = from->X.value(0);
  double y1 = from->Y.value(0);
  double z1 = from->Z.value(0);
  double x2 = to->X.value(0);
  double y2 = to->Y.value(0);
  double z2 = to->Z.value(0);

  x1 += from->r11*n1 + from->r12*e1 + from->r13*u1;
  y1 += from->r21*n1 + from->r22*e1 + from->r23*u1;
  z1 += from->r31*n1 + from->r32*e1 + from->r33*u1;

  x2 += to->r11*n2 + to->r12*e2 + to->r13*u2;
  y2 += to->r21*n2 + to->r22*e2 + to->r23*u2;
  z2 += to->r31*n2 + to->r32*e2 + to->r33*u2;

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


  if (from->N->active() && from->E->active() && from->U->active() && 
      to->  N->active() && to->  E->active() && to->  U->active()  )
    {
      Parameter** b = parlist.begin();

      *b++ = from->N;
      *b++ = from->E;
      *b++ = from->U;

      *b++ = to->N;
      *b++ = to->E;
      *b++ = to->U;
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
