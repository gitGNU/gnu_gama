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
 *  $Id: g3_model_distance.cpp,v 1.1 2004/01/25 11:07:13 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>


using namespace std;
using namespace GNU_gama::g3;


bool Model::revision_visit(Distance* d)
{
  if (!d->active()) return false;
  
  Point* from = points->find(d->from);
  Point* to   = points->find(d->to  );
  
  if ( from == 0       ||  to == 0        ) return d->set_active(false);
  if ( from->unused()  ||  to->unused()   ) return d->set_active(false);
  if (!from->has_blh() || !from->has_blh()) return d->set_active(false);

  active_obs->push_back(d);

  dm_rows++;                             // design matrix
  if (from->free_horizontal_position())
    {
      if (!from->N.index()) from->N_.set_index(++dm_cols);
      if (!from->E.index()) from->E_.set_index(++dm_cols);
      dm_floats += 2;
    }
  if (from->free_height())
    {
      if (!from->U.index()) from->U_.set_index(++dm_cols);
      dm_floats += 1;
    }

  if (to->free_horizontal_position())
    {
      if (!to->N.index()) to->N_.set_index(++dm_cols);
      if (!to->E.index()) to->E_.set_index(++dm_cols);
      dm_floats += 2;
    }
  if (to->free_height())
    {
      if (!to->U.index()) to->U_.set_index(++dm_cols);
      dm_floats += 1;
    }

  return d->active();
}


void Model::linearization_visit(Distance* d)
{
  const Point* from = points->find(d->from);
  const Point* to   = points->find(d->to  );

  // ....
}

