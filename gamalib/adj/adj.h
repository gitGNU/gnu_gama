/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: adj.h,v 1.1 2002/10/12 11:14:08 cepek Exp $
 */

#include <gamalib/exception.h>
#include <gamalib/sparse/smatrix.h>
#include <gamalib/sparse/sbdiagonal.h>
#include <gamalib/sparse/intlist.h>

#ifndef GaMaLib_Adj__adjustment_class__h
#define GaMaLib_Adj__adjustment_class__h


namespace GaMaLib {

  class AdjInputData {
  public:
    
    SparseMatrix <>  A;
    BlockDiagonal<>  cov;
    IntegerList  <>  minx;
    
  };
  


  class Adj {
    
    const AdjInputData* data;
    bool  solved;
    int   n_obs_, n_par_;
    
    void init(const AdjInputData*);

  public:
    
    enum algorithm {gso, svd};

    Adj () : data(0) { init(0); }
    ~Adj() { delete data; }

    int n_obs() const { return n_obs_; }
    int n_par() const { return n_par_; }

    void set(const AdjInputData* inp) { }
    void preferred_algorithm(Adj::algorithm);

  };
  
}  // namespace GaMaLib

#endif









