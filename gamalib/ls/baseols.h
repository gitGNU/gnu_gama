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
 *  $Id: baseols.h,v 1.2 2004/06/20 20:54:51 cepek Exp $
 */

#ifndef GaMa_BaseOLS_h
#define GaMa_BaseOLS_h

#include <gamalib/matvec.h>

namespace GaMaLib {


template <typename Float, typename Exc>
class BaseOLS {

public:
  BaseOLS() {}
  BaseOLS(const gMatVec::Mat<Float, Exc>& A, const gMatVec::Vec<Float, Exc>& b)
    : pA(&A), pb(&b), pw(0), is_solved(false) {}
  BaseOLS(const gMatVec::Mat<Float, Exc>& A, const gMatVec::Vec<Float, Exc>& b,
          const gMatVec::Vec<Float, Exc>& w)
    : pA(&A), pb(&b), pw(&w), is_solved(false) {}
  virtual ~BaseOLS() {}

  virtual void reset(const gMatVec::Mat<Float, Exc>& A, 
             const gMatVec::Vec<Float, Exc>& b) {
    pA = &A;
    pb = &b;
    pw = 0;
    is_solved = false;
  }
  virtual void reset(const gMatVec::Mat<Float, Exc>& A, 
             const gMatVec::Vec<Float, Exc>& b,
             const gMatVec::Vec<Float, Exc>& w)
  {
    pA = &A;
    pb = &b;
    pw = &w;
    is_solved = false;
  }

  const gMatVec::Vec<Float, Exc>& solve(gMatVec::Vec<Float, Exc>& x) 
    { 
      solve_me(); return x = BaseOLS::x; 
    }
  const gMatVec::Vec<Float, Exc>& solve() { solve_me(); return x; }
  const gMatVec::Vec<Float, Exc>& residuals(gMatVec::Vec<Float, Exc>& res) 
    { 
      solve_me(); return res = r; 
    }
  const gMatVec::Vec<Float, Exc>& residuals() { solve_me(); return r; }

  Float trwr();   // trans(r)*w*r
  virtual gMatVec::Index defect() = 0;


  virtual void  q_xx(gMatVec::Mat<Float, Exc>&);        // weight coefficients 
  virtual Float q_xx(gMatVec::Index, gMatVec::Index) = 0; // w. coeff. (xi,xj)
  virtual Float q_bb(gMatVec::Index, gMatVec::Index) = 0; //           (bi,bj)
  virtual Float q_bx(gMatVec::Index, gMatVec::Index) = 0; //           (bi,xj)

  virtual bool lindep(gMatVec::Index) = 0; // linearly dependent column/unknown
  virtual void min_x() = 0;
  virtual void min_x(gMatVec::Index, gMatVec::Index[]) = 0;

  virtual Float cond() { return 0; }  // condition number (0 if not available)

protected:

  // solve_me() must compute vectors x, r, sqrt_w and set is_solved=true
  virtual void solve_me() = 0;

  const gMatVec::Mat<Float, Exc>* pA;
  const gMatVec::Vec<Float, Exc>* pb;
  const gMatVec::Vec<Float, Exc>* pw;
  gMatVec::Vec<Float, Exc> x;
  gMatVec::Vec<Float, Exc> r;
  gMatVec::Vec<Float, Exc> sqrt_w;
  bool is_solved;

};

// ................................................................

template <typename Float, typename Exc>
Float BaseOLS<Float, Exc>::trwr()
{
  if (!is_solved) solve_me();

  Float s = 0, p;
  for (gMatVec::Index i = 1; i <= r.dim(); i++) {
      p = sqrt_w(i) * r(i);
      s += p*p;
  }

  return s;
}

template <typename Float, typename Exc>
void BaseOLS<Float, Exc>::q_xx(gMatVec::Mat<Float, Exc>& cxx)
{
  if (!is_solved) solve_me();

  const gMatVec::Index x_dim = x.dim();
  if (cxx.rows() != x_dim || cxx.cols() != x_dim)
    cxx.reset(x_dim, x_dim);

  for (gMatVec::Index i = 1; i <= x_dim; i++) {
    cxx(i,i) = q_xx(i, i);
    for (gMatVec::Index j = i+1; j <= x_dim; j++)
      cxx(i,j) = cxx(j,i) = q_xx(i, j);
  }
}

}
#endif

