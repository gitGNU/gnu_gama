/*  
    C++ Matrix/Vector templates (GNU GaMa / gMatVec 0.9.19)
    Copyright (C) 1999  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the gMatVec C++ Matrix/Vector template library.
    
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
 *  $Id: choldec.h,v 1.6 2002/07/11 20:54:15 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec_Cholesky_Decomposition__h_
#define gMatVec_Cholesky_Decomposition__h_

#include <gmatvec/vecbase.h>


namespace gMatVec {

/* Cholesky Decomposition of Positive Definite Matrix
 * ==================================================
 *
 * Two variants of implementation are possible:
 *
 *    a)   A = L*trans(L)
 *    b)   A = L*D*trans(L),   where D is diagonal and L has 1 on diagonal 
 */


template <class Float=double, class Exc=Exception>
class CholDec {

  Float  tol_;
  
public:

  CholDec(Float t=1e-8) : tol_(t) {}

  Float  cholTol() const  { return tol_; }
  Float  cholTol(Float t) { tol_ = t; return tol_; } 

  virtual void cholDec() = 0;    // `in situ' Cholesky decomposition
  virtual void solve(Vec<Float, Exc>& rhs) const = 0;
};


}   // namespace gMatVec

#endif



