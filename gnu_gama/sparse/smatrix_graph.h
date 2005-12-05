/*  
    Geodesy and Mapping C++ Library (GNU Gama)
    Copyright (C) 2004  Ales Cepek <cepek@gnu.org>

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
 *  $Id: smatrix_graph.h,v 1.1 2005/12/05 19:01:45 cepek Exp $
 */

#include <gnu_gama/sparse/smatrix.h>

#ifndef GNU_gama_matrix_graph_h___GNU_Gama_MatrixGraph
#define GNU_gama_matrix_graph_h___GNU_Gama_MatrixGraph


namespace GNU_gama {

  template <typename Float=double>
  class SparseMatrixGraph 
  {
  public:
    
    SparseMatrixGraph(const SparseMatrix<Float>* const m) : sparse(m) {}


    bool connected() const;

  private:
    
    const SparseMatrix<Float>* const sparse;

  };
}

#include <gnu_gama/sparse/smatrix_graph_connected.h>

#endif








