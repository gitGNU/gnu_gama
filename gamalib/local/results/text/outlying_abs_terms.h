/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 1999  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU GaMa / GaMaLib C++ Library.
    
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
 *  $Id: outlying_abs_terms.h,v 1.3 2003/03/01 18:16:18 cepek Exp $
 */

#ifndef GaMa_GaMaProg_Vybocujici_Absolutni_Cleny_h_
#define GaMa_GaMaProg_Vybocujici_Absolutni_Cleny_h_

#include <iomanip>
#include <gamalib/local/network.h>
#include <gamalib/local/pobs/format.h>
#include <gamalib/statan.h>

namespace GaMaLib {

template <class OutStream>
void OutlyingAbsoluteTerms(GaMaLib::LocalNetwork* IS, OutStream& out)
{
  using namespace std;
  using namespace GaMaLib;
  
  if (!IS->huge_abs_terms()) return;
  
  out << T_GaMa_abstrm_Review_of_outlying_abs_terms << "\n"
      << underline(T_GaMa_abstrm_Review_of_outlying_abs_terms, '*') << "\n\n";

  out.width(IS->maxw_obs());
  out << "i" << " ";
  out.width(IS->maxw_id());
  out << T_GaMa_standpoint << " ";
  out.width(IS->maxw_id());
  out << T_GaMa_target << T_GaMa_abstrm_header1;
  {  // for ...
    for (int i=0; i < (IS->maxw_obs() + 2*(IS->maxw_id()) + 13); i++) 
      out << "=";
  }  // for ...
  out << T_GaMa_abstrm_header2;
  out.flush();
  
  PointID predcs = "";   // previous standpoint ID
  for (int i=1; i<=IS->sum_observations(); i++)
    {
      Observation* pm = IS->ptr_obs(i);
      if (IS->test_abs_term(i))
        {
          out.width(IS->maxw_obs());
          out << i << " ";
          
          PointID cs = pm->from();
          out.width(IS->maxw_id());
          if (cs != predcs)
            out << cs.c_str();
          else
            out << " ";
          out << " ";
          
          PointID cc = pm->to();
          out.width(IS->maxw_id());
          out << cc.c_str();
          out.setf(ios::fixed, ios::floatfield);
          
          {   // ************************************************
            if (Distance* d = dynamic_cast<Distance*>(pm))
              {
                out << T_GaMa_distance;
                out.precision(5);
                out.width(12);
                Double m = d->value();
                out << m << " ";
              }
            else if (Direction* s = dynamic_cast<Direction*>(pm))
              {
                out << T_GaMa_direction;
                out.precision(6);
                out.width(12);
                Double m = R2G*(s->value());
                out << m << " ";
              }
            else if (Angle* u = dynamic_cast<Angle*>(pm))
              {
                out << '\n';
                out.width(IS->maxw_obs() + 2 + 2*(IS->maxw_id()));
                out << (u->fs()).c_str();
                out << T_GaMa_angle;
                out.precision(6);
                out.width(12);
                Double m = R2G*(u->value());
                out << m << " ";
              }
            else if (Z_Angle* z = dynamic_cast<Z_Angle*>(pm))
              {
                out << T_GaMa_z_angle;
                out.precision(6);
                out.width(12);
                Double m = R2G*(z->value());
                out << m << " ";
              }
            else if (S_Distance* p = dynamic_cast<S_Distance*>(pm))
              {
                out << T_GaMa_s_distance;
                out.precision(5);
                out.width(12);
                Double m = p->value();
                out << m << " ";
              }
            else if (H_Diff* h = dynamic_cast<H_Diff*>(pm))
              {
                out << T_GaMa_levell;
                out.precision(5);
                out.width(12);
                Double m = h->value();
                out << m << " ";
              }
            else
              {
                // vectors and coordinates have absoluted term 
                // always 0 in Network::test_abs_term()

                throw GaMaLib::Exception(
                   "GaMa internal error - unknown observation\n"
                   "$Id: outlying_abs_terms.h,v 1.3 2003/03/01 18:16:18 cepek Exp $");
              }
          }   // ************************************************
          
          out << setiosflags(ios::scientific) << setprecision(5);
          out << setw(13) << IS->rhs(i);       // 1.1.56 << pm->rhs();
          out << '\n';
          out.flush();
        }
    }
  
  out << "\n\n";
}

}

#endif

