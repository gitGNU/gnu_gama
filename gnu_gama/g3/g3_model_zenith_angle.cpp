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
 *  $Id: g3_model_zenith_angle.cpp,v 1.5 2004/08/18 13:00:57 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/radian.h>

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
  
  E_3 from_vertical, from_to, p1, v1, p2, v2;
  
  p1.set(from->X(), from->Y(), from->Z());
  v1 = from_vertical = vertical(from);  
  v1 *= z->from_dh;
  p1 += v1;                               // instrument
  
  p2.set(to->X(), to->Y(), to->Z());
  v2  = vertical(to);
  v2 *= z->to_dh;
  p2 += v2;                               // target
  
  from_to  = p2;
  from_to -= p1;                          // instrument --> target

  R_3 R;
  R.set_rotation(from->B(), from->L());   // dif_NEU --> dif_XYZ
  E_3 local;
  R.inverse(from_to, local);

  const double r = std::sqrt(local.e1*local.e1 + local.e2*local.e2);
  const double s = local.e1*local.e1 + local.e2*local.e2 + local.e3*local.e3;
  const double q = 1.0/(r*s);

  // pd - partial derivatives for the occupied station
  E_3 pd( -local.e1*q, -local.e2*q, r/s); 

  // nonzero derivatives in project equations
  A->new_row();
  if (from->free_horizontal_position())
    {
      A->add_element(pd.e1, from->N.index());
      A->add_element(pd.e2, from->E.index());
    }
  if (from->free_height())
    {
      A->add_element(pd.e3, from->U.index());
    }

  E_3 tmp;
  R.rotation(pd, tmp);
  tmp *= -1.0;
  R.set_rotation(to->B(),to->L());
  // pd - partial derivatives for the target station
  R.inverse (tmp, pd);   

  if (to->free_horizontal_position())
    {
      A->add_element(pd.e1, to->N.index());
      A->add_element(pd.e2, to->E.index());
    }
  if (to->free_height())
    {
      A->add_element(pd.e3, to->U.index());
    }


  // right hand site
  
  double za = angle(from_vertical, from_to); 

  /************************************************************/
  /*          !!! add refraction correction here !!!          */
  /************************************************************/
  std::cout << "??? zenith  " << (z->obs()*GON_TO_RAD - za)*RAD_TO_GON 
            << "\t" << z->obs() << "\t" << za*RAD_TO_GON 
            << "\n";
  rhs(++rhs_ind) = z->obs()*GON_TO_RAD - za;
}

