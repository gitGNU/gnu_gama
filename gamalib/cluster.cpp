/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: cluster.cpp,v 1.1 2001/12/07 12:22:42 cepek Exp $
 */

#include <gamalib/cluster.h>
#include <algorithm>

using namespace GaMaLib;


Cluster::~Cluster()
{
  for (ObservationList::iterator i=observation_list.begin(); 
                                 i!=observation_list.end() ;  ++i)
    {
      delete *i;
    }
}


void Cluster::update()
{
  act_count = 0;
  int index = 0;
  Observation* p;
  for (ObservationList::iterator i=observation_list.begin(); 
                                 i!=observation_list.end(); ++i)
    {
      p = (*i);
      p->cluster = this;
      p->cluster_index = index++;
      if (p->active()) act_count++;
    }
}

Cov Cluster::activeCov() const
{
  const Index M = covariance_matrix.rows();
  const Index B = covariance_matrix.bandWidth();
  const Index N = activeCount();
  Index temp = B;
  if (N-1 < B) temp = N-1;
  Cov C(N, temp);     // vc++ ... std::min<Index>(B, N-1)

  Index row = 1;
  Index col = row;
  for (Index i=0; i<M; ++i)
    if (observation_list[i]->active())
      {
        for (Index j=0; j<=B && i+j<M; ++j)
          {
            if (observation_list[i+j]->active())
              {
                C(row, col++) = covariance_matrix(i+1, i+j+1);
              }
          }
        col = ++row;
      }

  return C;
}





