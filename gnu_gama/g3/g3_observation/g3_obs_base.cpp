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
 *  $Id: g3_obs_base.cpp,v 1.8 2003/04/10 16:12:03 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation.h>


using namespace GNU_gama::g3;
using std::size_t;


void Observation::linearization(GNU_gama::SparseVector<>& row)
{
  prepare_to_linearization();
  
  row.reset();
  ParameterTree tree(parlist);
  ParameterTree::iterator b=tree.begin(), e=tree.end();

  double     d; 
  Parameter* p;

  while (b != e)
    {
      p = *b;
      if (p->free())
        {
          d = derivative(p);
          if (d) row.add(p->index(), d);
        }
      ++b;
    }
}


double Observation::numerical_derivative(Parameter* p) 
{
  if (!p->free()) return 0;
 
  double p_correction = p->correction();
  double d1, d2, d, h;
  { 
    // L'4(x) = +2/24*y(x-2h)  -4/6*y(x-h)  +4/6*y(x+h) -2/24*y(x+2h)

    d = p->step_size() + p_correction;
    h = d - p_correction;               // temp = x+h; h = temp-x 

    p->set_correction(p_correction - 2*h);  d2  = parlist_value();
    p->set_correction(p_correction + 2*h);  d2 -= parlist_value();  
    p->set_correction(p_correction +  h );  d1  = parlist_value();
    p->set_correction(p_correction -  h );  d1 -= parlist_value();

    d = (d2 + 8.0*d1) / (12.0*h);
  }
  p->set_correction(p_correction);

  return d;
}
