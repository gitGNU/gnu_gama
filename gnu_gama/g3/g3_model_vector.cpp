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
 *  $Id: g3_model_vector.cpp,v 1.5 2004/03/24 19:27:07 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>

using namespace GNU_gama::g3;

bool Model::revision_visit(Vector* v)
{
  if (!v->active()) return false;

  Point* from = points->find(v->from);
  Point* to   = points->find(v->to  );

  if ( from == 0      ||  to == 0      ) return v->set_active(false);
  if ( from->unused() ||  to->unused() ) return v->set_active(false);
  if (!from->has_position()            ) return v->set_active(false);
  if (!to  ->has_position()            ) return v->set_active(false);

  active_obs->push_back(v);

  update_index(from->N_);
  update_index(from->E_);
  update_index(from->U_);
  update_index(to  ->N_);
  update_index(to  ->E_);
  update_index(to  ->U_);

  dm_rows += v->dimension();            // design matrix

  if (from->free_horizontal_position())  dm_floats += 6;
  if (from->free_height())               dm_floats += 3;
  if (to->free_horizontal_position())    dm_floats += 6;
  if (to->free_height())                 dm_floats += 3;

  return v->active();
}


void Model::linearization_visit(Vector* v)
{
  Point* from = points->find(v->from);
  Point* to   = points->find(v->to  );
 
  for (int i=1; i<=3; i++)
   {
     double tx=0, ty=0, tz=0;
     switch (i) {
     case 1: tx = 1.0; break;
     case 2: ty = 1.0; break;
     case 3: tz = 1.0; break;
     };

     from->set_diff_XYZ(-tx, -ty, -tz);
     to  ->set_diff_XYZ( tx,  ty,  tz);
  
     // nonzero derivatives in project equations
     A->new_row();

     if (from->free_horizontal_position())
       {
         A->add_element(from->diff_N(), from->N.index());
         A->add_element(from->diff_E(), from->E.index());
       }
     if (from->free_height())
       {
         A->add_element(from->diff_U(), from->U.index());
       }
     
     if (to->free_horizontal_position())
       {
         A->add_element(to->diff_N(), to->N.index());
         A->add_element(to->diff_E(), to->E.index());
       }
     if (to->free_height())
       {
         A->add_element(to->diff_U(), to->U.index());
       }
   }

   // right hand site
   {
     double dx = to->X_dh(v->to_dh) - from->X_dh(v->from_dh);
     double dy = to->Y_dh(v->to_dh) - from->Y_dh(v->from_dh);
     double dz = to->Z_dh(v->to_dh) - from->Z_dh(v->from_dh);
 
     rhs(++rhs_ind) = dx - v->dx();
     rhs(++rhs_ind) = dy - v->dy();
     rhs(++rhs_ind) = dz - v->dz();
   }
}

