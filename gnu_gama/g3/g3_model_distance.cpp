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
 *  $Id: g3_model_distance.cpp,v 1.2 2004/01/26 19:03:09 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>

using namespace GNU_gama::g3;

//---  namespace {
//---  
//---    struct distance
//---    {
//---      Point   *A, *B;
//---      double   Na, Ea, Ua,  Nb, Eb, Ub;
//---  
//---      distance(Point* a, Point* b, double ua, double ub)
//---        : A(a), B(b), Na(0.), Ea(0.), Ua(ua), Nb(0.), Eb(0.), Ub(ub)
//---      {
//---      }
//---  
//---      double value()
//---      {
//---        double dx 
//---          = B->X() + B->x_transform(Nb, Eb, Ub)
//---          - A->X() - A->x_transform(Na, Ea, Ua);
//---  
//---        double dy 
//---          = B->Y() + B->y_transform(Nb, Eb, Ub)
//---          - A->Y() - A->y_transform(Na, Ea, Ua);
//---  
//---        double dz 
//---          = B->Z() + B->z_transform(Nb, Eb, Ub)
//---          - A->Z() - A->z_transform(Na, Ea, Ua);
//---  
//---        return std::sqrt(dx*dx + dy*dy + dz*dz);
//---      }
//---  
//---      double d(double& p)
//---      {
//---        const double diff = 0.5;
//---        double dd, tmp = p;
//---  
//---        p  += diff;
//---        dd  = value();
//---        p  -= 2*diff;
//---        dd -= value();
//---        dd /= (2*diff);
//---  
//---        p = tmp;
//---        return dd;
//---      }
//---  
//---    };
//---  }


bool Model::revision_visit(Distance* d)
{
  if (!d->active()) return false;
  
  Point* from = points->find(d->from);
  Point* to   = points->find(d->to  );
  
  if ( from == 0      ||  to == 0      ) return d->set_active(false);
  if ( from->unused() ||  to->unused() ) return d->set_active(false);
  if (!from->has_position()            ) return d->set_active(false);
  if (!to  ->has_position()            ) return d->set_active(false);

  active_obs->push_back(d);

  dm_rows += d->dimension();            // design matrix
  if (from->free_horizontal_position())
    {
      update_index(from->N_);
      update_index(from->E_);
      dm_floats += 2;
    }
  if (from->free_height())
    {
      update_index(from->U_);
      dm_floats += 1;
    }

  if (to->free_horizontal_position())
    {
      update_index(to->N_);
      update_index(to->E_);
      dm_floats += 2;
    }
  if (to->free_height())
    {
      update_index(to->U_);
      dm_floats += 1;
    }

  return d->active();
}


void Model::linearization_visit(Distance* d)
{
  Point* from = points->find(d->from);
  Point* to   = points->find(d->to  );


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



  // nonzero derivatives in proeject equations
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



  // right hand site
  {
    double dx = to->X_dh(d->to_dh) - from->X_dh(d->from_dh);
    double dy = to->Y_dh(d->to_dh) - from->Y_dh(d->from_dh);
    double dz = to->Z_dh(d->to_dh) - from->Z_dh(d->from_dh);
    double D  = std::sqrt(dx*dx + dy*dy + dz*dz);

    rhs(++rhs_ind) = D - d->obs();
  }
}

