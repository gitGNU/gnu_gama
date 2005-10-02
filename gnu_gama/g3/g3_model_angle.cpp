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
 *  $Id: g3_model_angle.cpp,v 1.2 2005/10/02 16:00:27 cepek Exp $
 */

#include <gnu_gama/g3/g3_model.h>
#include <gnu_gama/radian.h>

using namespace GNU_gama::g3;

bool Model::revision_visit(Angle* angle)
{
  if (!angle->active()) return false;

  Point* from  = points->find(angle->from);
  Point* left  = points->find(angle->left);
  Point* right = points->find(angle->right);

  if ( from  == 0      ) return angle->set_active(false);
  if ( from->unused()  ) return angle->set_active(false);
  if (!from->has_xyz() ) return angle->set_active(false);

  if ( left  == 0      ) return angle->set_active(false);
  if ( left->unused()  ) return angle->set_active(false);
  if (!left->has_xyz() ) return angle->set_active(false);

  if ( right  == 0     ) return angle->set_active(false);
  if ( right->unused() ) return angle->set_active(false);
  if (!right->has_xyz()) return angle->set_active(false);

  active_obs->push_back(angle);

  update_index(from->N);
  update_index(from->E);
  update_index(from->U);
  update_index(left->N);
  update_index(left->E);
  update_index(left->U);
  update_index(right->N);
  update_index(right->E);
  update_index(right->U);

  dm_rows += angle->dimension();      // design matrix
   
  if (from ->free_position()) dm_floats += 2;
  if (from ->free_height  ()) dm_floats += 1;
  if (left ->free_position()) dm_floats += 2;
  if (left ->free_height  ()) dm_floats += 1;
  if (right->free_position()) dm_floats += 2;
  if (right->free_height  ()) dm_floats += 1;

  return angle->active();
}


void Model::linearization_visit(Angle* pangle)
{
  Point* from  = points->find(pangle->from);
  Point* left  = points->find(pangle->left);
  Point* right = points->find(pangle->right);
  
  const double sca = Angular().scale();
  const double scl = Linear ().scale();
  const double sc  = sca / scl;

  E_3 FI = instrument(from,  pangle->from_dh);
  E_3 FL = instrument(left,  pangle->left_dh); 
  E_3 FR = instrument(right, pangle->right_dh); 
  E_3 FV = vertical(from);
  FL    -= FI;
  FR    -= FI;
  double ql  = sqrt(FL.e1*FL.e1 + FL.e2*FL.e2 + FL.e3*FL.e3);
  double qr  = sqrt(FR.e1*FR.e1 + FR.e2*FR.e2 + FR.e3*FR.e3);
  if (ql) ql = 1.0/ql;
  if (qr) qr = 1.0/qr;
  FL *= ql;
  FR *= qr;
  E_3 VL, VR;
  VL.cross(FV, FL);
  VR.cross(FV, FR);


  R_3 tran;
  tran.set_rotation(from->B(), from->L());

  E_3 Lxyz(left ->X.init_value() - from->X.init_value(), 
           left ->Y.init_value() - from->Y.init_value(), 
           left ->Z.init_value() - from->Z.init_value());
  E_3 Rxyz(right->X.init_value() - from->X.init_value(), 
           right->Y.init_value() - from->Y.init_value(), 
           right->Z.init_value() - from->Z.init_value());

  E_3 Lneu, Rneu;
  tran.inverse(Lxyz, Lneu);
  tran.inverse(Rxyz, Rneu);

  const double dl = sqrt (Lneu.e1*Lneu.e1 + Lneu.e2*Lneu.e2);
  const double dr = sqrt (Rneu.e1*Rneu.e1 + Rneu.e2*Rneu.e2);
  const double sl = atan2(Lneu.e2, Lneu.e1);
  const double sr = atan2(Rneu.e2, Rneu.e1);

  const double psl = sin(sl)/dl;
  const double pcl = cos(sl)/dl;
  const double psr = sin(sr)/dr;
  const double pcr = cos(sr)/dr;

  E_3 Lcoef(-psl,  pcl, 0.0);
  E_3 Rcoef(-psr,  pcr, 0.0);

  E_3 Fcoef( Rcoef );  Fcoef -= Lcoef;  Fcoef *= -1.0;

  tran.rotation(Lcoef, Lxyz);
  tran.rotation(Rcoef, Rxyz);
  tran.set_rotation(left->B(), left->L());
  tran.inverse(Lxyz, Lcoef);
  tran.set_rotation(right->B(), right->L());
  tran.inverse(Rxyz, Rcoef);


  // nonzero derivatives in project equations
  A->new_row();

  if (from->N.free()) A->add_element(Fcoef.e1*sc, from->N.index());
  if (from->E.free()) A->add_element(Fcoef.e2*sc, from->E.index());
  if (from->U.free()) A->add_element(Fcoef.e3*sc, from->U.index());

  if (left->N.free()) A->add_element(Lcoef.e1*sc, left->N.index());
  if (left->E.free()) A->add_element(Lcoef.e2*sc, left->E.index());
  if (left->U.free()) A->add_element(Lcoef.e3*sc, left->U.index());

  if (right->N.free()) A->add_element(Rcoef.e1*sc, right->N.index());
  if (right->E.free()) A->add_element(Rcoef.e2*sc, right->E.index());
  if (right->U.free()) A->add_element(Rcoef.e3*sc, right->U.index());


  // right hand site

  const double angle = GNU_gama::angle(VL, VR);
  rhs(++rhs_ind) = (pangle->obs() - angle)*sca;
}

