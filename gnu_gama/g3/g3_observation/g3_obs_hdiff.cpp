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
 *  $Id: g3_obs_hdiff.cpp,v 1.7 2003/04/08 16:41:51 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation.h>
#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/g3/g3_observation/g3_der_hdiff.h>

using namespace GNU_gama::g3;


double HeightDiff::parlist_value() const
{
  Parameter** p = parlist.begin();
  double s = (*p++)->value(time);
  double t = (*p++)->value(time);

  return t - s;
}


void HeightDiff::parlist_init(Model* model)
{
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!! this is not a real solution. used here just for testing !!! 
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ellipsoid = &model->ellipsoid;
  
  if (!active())  return;
  
  Point* from = model->points.find(name[0]);
  Point* to   = model->points.find(name[1]);
  
  if (from == 0 || to == 0)
    {
      set_active(false);
      return;
    }
  
  
  if (from->U->active() && to->U->active())
    {
      Parameter** b = parlist.begin();
      *b++ = from->U;
      *b++ = to->U;
    }
  else
    {
      set_active(false);
      return;
    }
}


double HeightDiff::derivative(Parameter* p)
{
  HeightDiffAnalyticalDerivative* ad = 
    dynamic_cast<HeightDiffAnalyticalDerivative*>(p);
  
  if (ad)
    return ad->analytical_derivative(this);
  else
    return numerical_derivative(p);
}
