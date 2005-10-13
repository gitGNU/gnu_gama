/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2005  Ales Cepek <cepek@gnu.org>

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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 *  $Id: g3_model_height.cpp,v 1.4 2005/10/13 18:57:50 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/radian.h>

using namespace GNU_gama::g3;

bool Model::revision_visit(Height* height)
{
  if (!height->active()) return false;
  
  Point* point = points->find(height->id);
  
  if ( point == 0                ) return height->set_active(false);
  if ( point->unused()           ) return height->set_active(false);
  if (!point->test_model_height()) return height->set_active(false);

  active_obs->push_back(height);

  update_index(point->U);
  
  dm_rows += height->dimension();      // design matrix

  if (point->free_height()) dm_floats += 1;

  return height->active();
}


void Model::linearization_visit(Height* height)
{
  Point* point = points->find(height->id);
  
  // nonzero derivatives in project equations
  A->new_row();

  if (point->free_height())  A->add_element(1, point->U.index());

  // right hand site

  rhs(++rhs_ind) = (height->obs() - point->model_height())*Linear().scale();
}

