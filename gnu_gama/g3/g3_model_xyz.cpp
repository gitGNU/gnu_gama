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
 *  $Id: g3_model_xyz.cpp,v 1.4 2004/04/23 22:01:31 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>

using namespace GNU_gama::g3;

bool Model::revision_visit(XYZ* xyz)
{
  if (!xyz->active()) return false;

  Point* point = points->find(xyz->id);

  if ( point == 0           ) return xyz->set_active(false);
  if ( point->unused()      ) return xyz->set_active(false);
  if (!point->has_position()) return xyz->set_active(false);

  active_obs->push_back(xyz);

  update_index(point->N);
  update_index(point->E);
  update_index(point->U);

  dm_rows += xyz->dimension();            // design matrix
  if (point->free_horizontal_position())   dm_floats += 6;
  if (point->free_height())                dm_floats += 3;

  return xyz->active();
}


void Model::linearization_visit(XYZ* xyz)
{
  Point* point = points->find(xyz->id);
 
  for (int i=1; i<=3; i++)
   {
     double tx=0, ty=0, tz=0;
     switch (i) {
     case 1: tx = 1.0; break;
     case 2: ty = 1.0; break;
     case 3: tz = 1.0; break;
     };

     point->set_diff_XYZ(tx, ty, tz);
  
     // nonzero derivatives in project equations
     A->new_row();

     if (point->free_horizontal_position())
       {
         A->add_element(point->diff_N(), point->N.index());
         A->add_element(point->diff_E(), point->E.index());
       }
     if (point->free_height())
       {
         A->add_element(point->diff_U(), point->U.index());
       }     
   }

   // right hand site
   {
     double x = point->X();
     double y = point->Y();
     double z = point->Z();
   
     rhs(++rhs_ind) = xyz->x() - x;
     rhs(++rhs_ind) = xyz->y() - y;
     rhs(++rhs_ind) = xyz->z() - z;
   }
}

