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
 *  $Id: g3_obs_hdiff.cpp,v 1.2 2003/03/23 18:39:53 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation.h>
#include <gnu_gama/g3/g3_model.h>

using namespace GNU_gama::g3;


void H_diff::init_parameters(Model* model)
{
  if (!active())  return;
 
  Point* from = model->points.find(name[0]);
  Point* to   = model->points.find(name[1]);

  if (from == 0 || to == 0)
    {
      set_active(false);
      return;
    }

  if (from->state(Point::height) && to->state(Point::height))
    {
      Parameter** b = parlist.begin();
      *b++ = from->H;
      *b++ = to->H;
    }
  else
    {
      set_active(false);
      return;
    }
}

