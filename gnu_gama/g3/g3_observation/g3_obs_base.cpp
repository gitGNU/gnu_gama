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
 *  $Id: g3_obs_base.cpp,v 1.3 2003/03/26 17:33:47 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation.h>


using namespace GNU_gama::g3;
using std::size_t;


void Observation::linearization(GNU_gama::SparseVector<>& row)
{
  row.reset();
  ParameterTree tree(parlist);
  ParameterTree::iterator b=tree.begin(), e=tree.end();
  Parameter* p;
  size_t     n;
  double     d;
 
  while (b != e)
    {
      p = *b;
      if ((d = numerical_derivative(p)))
      {
        static int k=0;
        k++;
        row.add(k, d);
      }
      ++b;
    }
}


double Observation::numerical_derivative(Parameter* p) 
{
  if (p->fixed()) return 0;
 
  static double x = 1;
  x += 0.01;
  return x;
}
