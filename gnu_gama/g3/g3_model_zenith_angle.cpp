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
 *  $Id: g3_model_zenith_angle.cpp,v 1.1 2004/05/12 18:29:22 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>

using namespace GNU_gama::g3;

bool Model::revision_visit(ZenithAngle* z)
{
  if (!z->active()) return false;
  
  Point* from = points->find(z->from);
  Point* to   = points->find(z->to  );
  
  if ( from == 0      ||  to == 0      ) return z->set_active(false);
  if ( from->unused() ||  to->unused() ) return z->set_active(false);
  if (!from->has_position()            ) return z->set_active(false);
  if (!to  ->has_position()            ) return z->set_active(false);

  active_obs->push_back(z);

  update_index(from->N);
  update_index(from->E);
  update_index(from->U);
  update_index(to  ->N);
  update_index(to  ->E);
  update_index(to  ->U);
  
  dm_rows += z->dimension();            // design matrix

  if (from->free_horizontal_position())  dm_floats += 2;
  if (from->free_height())               dm_floats += 1;
  if (to  ->free_horizontal_position())  dm_floats += 2;
  if (to  ->free_height())               dm_floats += 1;

  return z->active();
}


void Model::linearization_visit(ZenithAngle* z)
{
  Point* from = points->find(z->from);
  Point* to   = points->find(z->to  );

  /* ++++  DEBUG  ++++++++++++++++++++++++++++++++++++++ */
  {
    double dx = to->X() - from->X();
    double dy = to->Y() - from->Y();
    double dz = to->Z() - from->Z();
    if (double dd = std::sqrt(dx*dx + dy*dy + dz*dz))
      {
        dx /= dd;
        dy /= dd;
        dz /= dd;
      }
    from->set_diff_XYZ(-dx, -dy, -dz);
    to  ->set_diff_XYZ( dx,  dy,  dz);
  }
  /* ----  DEBUG  --------------------------------------- */



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



  // right hand site  ... DEBUG !!!
  {
    double dx = to->X_dh(z->to_dh) - from->X_dh(z->from_dh);
    double dy = to->Y_dh(z->to_dh) - from->Y_dh(z->from_dh);
    double dz = to->Z_dh(z->to_dh) - from->Z_dh(z->from_dh);
    double D  = std::sqrt(dx*dx + dy*dy + dz*dz);

    rhs(++rhs_ind) = z->obs() - D;
  }
}

