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
 *  $Id: g3_model_height_diff.cpp,v 1.2 2005/07/27 15:11:08 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/radian.h>

using namespace GNU_gama::g3;

bool Model::revision_visit(HeightDiff* dh)
{
  if (!dh->active()) return false;
  
  Point* from = points->find(dh->from);
  Point* to   = points->find(dh->to  );
  
  if ( from == 0      ||  to == 0      ) return dh->set_active(false);
  if ( from->unused() ||  to->unused() ) return dh->set_active(false);
  if (!from->has_height()              ) return dh->set_active(false);
  if (!to  ->has_height()              ) return dh->set_active(false);

  active_obs->push_back(dh);

  update_index(from->N);
  update_index(from->E);
  update_index(from->U);
  update_index(to  ->N);
  update_index(to  ->E);
  update_index(to  ->U);
  
  dm_rows += dh->dimension();            // design matrix

  if (from->free_height())               dm_floats += 1;
  if (to  ->free_height())               dm_floats += 1;

  return dh->active();
}


void Model::linearization_visit(HeightDiff* dh)
{
  Point* from = points->find(dh->from);
  Point* to   = points->find(dh->to  );
  
  // nonzero derivatives in project equations
  A->new_row();

  if (from->free_height())  A->add_element(-1, from->U.index());
  if (to  ->free_height())  A->add_element(+1, to  ->U.index());

  // right hand site
 
  double h = to->height() - from->height();  
  // *** add here : correction for vertical deflections

  rhs(++rhs_ind) = dh->obs() - h;
}

