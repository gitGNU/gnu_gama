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
 *  $Id: g3_cluster_vec.cpp,v 1.4 2003/06/08 08:11:13 cepek Exp $
 */

#include <gnu_gama/g3/g3_observation/g3_cluster_vec.h>
#include <iomanip>

using namespace GNU_gama::g3;
using namespace std;


Vectors::~Vectors()
{
  for (GNU_gama::List<Vector*>::iterator 
         b=vectors.begin(), e=vectors.end(); b != e;  ++b)
    {
      delete *b;
    }
}

void Vectors::add(Vector* v)
{
  vectors.push_back(v);

  observation_list.push_back( new DiffX(v) );
  observation_list.push_back( new DiffY(v) );
  observation_list.push_back( new DiffZ(v) );
}


void Vectors::write_xml(std::ostream& out) const
{
  const bool single = vectors.size() == 1;

  if (!single)
    {
      out << "\n\t< multiple vectors cluster not implemented yet />\n\n";
      return;
    }
    
  const Vector* v = *vectors.begin();

  for (GNU_gama::List<Vector*>::const_iterator
       b = vectors.begin(), e=vectors.end();  b!=e;  ++b)
    {
      out.precision(5);
      out.setf(ios::fixed, ios::floatfield);

      out << "\n<vector>\n\t"
          << "<from>" << v->name[0] << "</from> "
          << "<to>"   << v->name[1] << "</to>\n\t"
          << "<dx>"   << v->dx() << "</dx> "
          << "<dy>"   << v->dy() << "</dy> "
          << "<dz>"   << v->dz() << "</dz>\n";

      if (single) 
        {
          out.precision(6);
          out.setf(ios::scientific, ios::floatfield);
     
          out << "\t"
              << "<cxx>" << covariance_matrix(1,1) << "</cxx> "
              << "<cxy>" << covariance_matrix(1,2) << "</cxy> "
              << "<cxz>" << covariance_matrix(1,3) << "</cxz>\n\t"
              << "<cyy>" << covariance_matrix(2,2) << "</cyy> "
              << "<cyz>" << covariance_matrix(2,3) << "</cyz>\n\t"
              << "<czz>" << covariance_matrix(3,3) << "</czz>\n";
        }

      out << "</vector>\n";
    }
}


void Vectors::parlist_init(Model* model)
{
  for (GNU_gama::List<Vector*>::iterator
       i = vectors.begin(), e = vectors.end();  i != e;  ++i)
    {
      (*i)->parlist_init(model);
    }
  
}
