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
 *  $Id: olssvd.h,v 1.1 2001/12/07 12:56:38 cepek Exp $
 */

#ifndef GaMa_OLS_svd_h
#define GaMa_OLS_svd_h

#include <gamalib/ls/baseols.h>
#include <cmath>

namespace GaMaLib {
  
template <class Float, class Exc>
class OLSsvd : virtual public BaseOLS<Float, Exc> {

  gMatVec::SVD<Float, Exc> svd;

public:
  OLSsvd() {}
  OLSsvd(const gMatVec::Mat<Float, Exc>& A, const gMatVec::Vec<Float, Exc>& b)
    : BaseOLS<Float, Exc>(A, b) {}
  OLSsvd(const gMatVec::Mat<Float, Exc>& A, const gMatVec::Vec<Float, Exc>& b,
         const gMatVec::Vec<Float, Exc>& w) : BaseOLS<Float, Exc>(A, b, w) {}
  
  void reset(const gMatVec::Mat<Float, Exc>& A, 
             const gMatVec::Vec<Float, Exc>& b)
    {
      BaseOLS<Float, Exc>::reset(A, b);
      svd.reset(A);
    }
  void reset(const gMatVec::Mat<Float, Exc>& A, 
             const gMatVec::Vec<Float, Exc>& b,
             const gMatVec::Vec<Float, Exc>& w)
    {
      BaseOLS<Float, Exc>::reset(A, b, w);
      svd.reset(A);
    }
  
  const gMatVec::Vec<Float, Exc>& solve(gMatVec::Vec<Float, Exc>& x)
    {
      return x = BaseOLS<Float, Exc>::solve();
    }
  const gMatVec::Vec<Float, Exc>& solve() { return BaseOLS<Float, Exc>::solve(); }
  
  gMatVec::Index defect() { return svd.nullity(); }
  bool  lindep(gMatVec::Index i) { return svd.lindep(i); }
  
  void  q_xx(gMatVec::Mat<Float, Exc>& C) { BaseOLS<Float, Exc>::q_xx(C); }
  Float q_xx(gMatVec::Index i, gMatVec::Index j)
    {
      if(!is_solved) solve_me();
      return svd.q_xx(i, j);
    }
  Float q_bb(gMatVec::Index i, gMatVec::Index j)
    {
      if (!is_solved) solve_me();
      return svd.q_bb(i, j) / (sqrt_w(i) * sqrt_w(j));
    }
  Float q_bx(gMatVec::Index i, gMatVec::Index j)
    {
      if (!is_solved) solve_me();
      return svd.q_bx(i, j) / sqrt_w(i);
    }
  
  void min_x()   {  svd.min_x(); }
  void min_x(gMatVec::Index n, gMatVec::Index x[]) { svd.min_x(n, x); }

  Float cond();
  
protected:
  
   void solve_me();
   
};

// ...................................................................

template <class Float, class Exc>
void OLSsvd<Float, Exc>::solve_me()
{
   if (is_solved) return;
   using namespace gMatVec; 

   sqrt_w.reset(pb->dim());
   if (pw)
   { 
      for (gMatVec::Index n = 1; n <= sqrt_w.dim(); n++)
         sqrt_w(n) = sqrt((*pw)(n));
      svd.reset(*pA, sqrt_w);                 // protected ==> public
      gMatVec::Vec<Float, Exc> pbw = *pb;
      for (gMatVec::Index i = 1; i <= pb->dim(); i++)
         pbw(i) *= sqrt_w(i);
      svd.solve(pbw, x);
   }
   else
   {
      svd.reset(*pA);
      for (gMatVec::Index n = 1; n <= sqrt_w.dim(); n++)
         sqrt_w(n) = 1;
      svd.solve(*pb, x);
   }
   r = *pA * x;
   r -= *pb;

   is_solved = true;
}

template <class Float, class Exc>
Float OLSsvd<Float, Exc>::cond()
{
  const Vec& W = svd.SVD_W();

  Float  f, sv_min=W(1), sv_max=W(1);

  for (Index i=2; i<=W.dim(); i++)
    if (!svd.lindep(i))
      {
        f = W(i); if (f < 0) f = -f;

        if (f < sv_min) sv_min = f;
        if (f > sv_max) sv_max = f;
      }

  return sv_max/sv_min;
}

}
#endif









