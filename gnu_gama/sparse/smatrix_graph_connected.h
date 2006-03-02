/*  
    Geodesy and Mapping C++ Library (GNU Gama)
    Copyright (C) 2005  Ales Cepek <cepek@gnu.org>

    This file is part of the GNU Gama C++ Library.
    
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
 *  $Id: smatrix_graph_connected.h,v 1.3 2006/03/02 13:41:21 cepek Exp $
 */


#include <matvec/matvec.h>

template <typename Float>
bool GNU_gama::SparseMatrixGraph<Float>::connected() const
{
  const Index cols = sparse->columns();
  const Index rows = sparse->rows();

  if (cols == 0 || rows == 0) return true;

  Vec<int> col(cols);
  Vec<int> row(rows);

  col.set_zero();
  for (Index i=1; i<=rows; i++) row(i) = i;

  col(1) = 1;     // there is always at least one unknown parameter

  bool  updated;
  Index r, row_count = rows;
  do {
    updated = false;
    
    r = 0;
    while (++r <= row_count)
      {
        const Index* s = sparse->ibegin(row(r));
        const Index* e = sparse->iend(row(r));

        for (const Index* i=s; i!=e; i++)
          if (col(*i))
          {
            for (const Index* i=s; i!=e; i++)  col(*i) = 1;
            
            updated = true;
            row(r--) = row(row_count--);
            break;
          }
      } 

  } while (updated);


  return row_count == 0;
}
