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
 *  $Id: coordinates.h,v 1.1 2001/12/07 13:02:30 cepek Exp $
 */

#ifndef GaMaLib_XML_adjusted_coordinates____h
#define GaMaLib_XML_adjusted_coordinates____h

#include <gamalib/local/gamadata.h>
#include <gamalib/local/network.h>
#include <gamalib/cluster.h>

template <class OutStream> void 
XML_adjusted_coordinates(GaMaLib::LocalNetwork* netinfo, 
                         OutStream& out, bool covm)
{
  using namespace std;
  using namespace GaMaLib;
  // using GaMaLib::Double;
  
  {
    bool sour = false;
    for (int i=1; i<=netinfo->sum_unknowns(); i++)
      if (netinfo->unknown_type(i) == 'X' || netinfo->unknown_type(i) == 'Z')
        {
          sour = true;
          break;
       }
    if (!sour) return;
  }

  
  out << "\n<coordinates>\n\n";

  // points ....
 
  const int maxdim = netinfo->sum_unknowns() + 1;
  Index* ind = new Index[maxdim];
  Index dim = 0;
  {
    const Vec& x = netinfo->solve();
    out.setf(ios::fixed, ios::floatfield);
    out.precision(5);
    for (PointData::const_iterator 
           i=netinfo->PD.begin(); i!=netinfo->PD.end(); ++i)
      {
        const Point& b = (*i).second;
        if (!b.active_xy()) continue;
        if (b.index_x()) 
          {
            ind[++dim] = b.index_x();
            ind[++dim] = b.index_y();
          }
        if (b.index_z()) 
          {
            ind[++dim] = b.index_z();
          }

        out << "<point id=\"" << (*i).first << "\"";
        if (b.index_x()) 
          {
            out << " x=\"" << (b.x()+x(b.index_x())/1000) << "\""
                << " y=\"" << (b.y()+x(b.index_y())/1000) << "\"";
          }        
        if (b.index_z()) 
          {
            out << " z=\"" << (b.z()+x(b.index_z())/1000) << "\"";
          }        
        out << " />\n";
      }
  }

  // covariances ....

  if (covm) 
    {
      out << "\n<cov-mat dim=\"" << dim << "\" band=\"" << dim-1 << "\" >\n\n";

      out.setf(ios::scientific, ios::floatfield);
      out.precision(7);
      const double m2 = netinfo->m_0() * netinfo->m_0();
      for (Index k=0, i=1; i<=dim; i++)
        for (Index j=i; j<=dim; j++)
          {
            out << m2*netinfo->qxx(ind[i], ind[j]);
            if (++k == 5)
              {
                k = 0;
                out << "\n";
              }
            else
              {
                out << " ";
              }
          }
      
      out << "\n\n</cov-mat>\n";
    }
  delete[] ind;

  out << "\n</coordinates>\n";  
}

#endif
