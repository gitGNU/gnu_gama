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
 *  $Id: olssvd.h,v 1.6 2005/03/27 17:43:26 cepek Exp $
 */

#ifndef GaMa_OLS_svd_h
#define GaMa_OLS_svd_h

#include <gamalib/ls/baseols.h>
#include <cmath>

namespace GaMaLib {
  
template <typename Float, typename Exc>
class OLSsvd : virtual public BaseOLS<Float, Exc> {

  GNU_gama::SVD<Float, Exc> svd;

public:
  OLSsvd() {}
  OLSsvd(const GNU_gama::Mat<Float, Exc>& A, const GNU_gama::Vec<Float, Exc>& b)
    : BaseOLS<Float, Exc>(A, b) {}
  OLSsvd(const GNU_gama::Mat<Float, Exc>& A, const GNU_gama::Vec<Float, Exc>& b,
         const GNU_gama::Vec<Float, Exc>& w) : BaseOLS<Float, Exc>(A, b, w) {}
  
  void reset(const GNU_gama::Mat<Float, Exc>& A, 
             const GNU_gama::Vec<Float, Exc>& b)
    {
      BaseOLS<Float, Exc>::reset(A, b);
      svd.reset(A);
    }
  void reset(const GNU_gama::Mat<Float, Exc>& A, 
             const GNU_gama::Vec<Float, Exc>& b,
             const GNU_gama::Vec<Float, Exc>& w)
    {
      BaseOLS<Float, Exc>::reset(A, b, w);
      svd.reset(A);
    }
  
  const GNU_gama::Vec<Float, Exc>& solve(GNU_gama::Vec<Float, Exc>& x)
    {
      return x = BaseOLS<Float, Exc>::solve();
    }
  const GNU_gama::Vec<Float, Exc>& solve() { return BaseOLS<Float, Exc>::solve(); }
  
  GNU_gama::Index defect() { return svd.nullity(); }
  bool  lindep(GNU_gama::Index i) { return svd.lindep(i); }
  
  void  q_xx(GNU_gama::Mat<Float, Exc>& C) { BaseOLS<Float, Exc>::q_xx(C); }
  Float q_xx(GNU_gama::Index i, GNU_gama::Index j)
    {
      if(!this->is_solved) solve_me();
      return svd.q_xx(i, j);
    }
  Float q_bb(GNU_gama::Index i, GNU_gama::Index j)
    {
      if (!this->is_solved) solve_me();
      return svd.q_bb(i, j) / (this->sqrt_w(i) * this->sqrt_w(j));
    }
  Float q_bx(GNU_gama::Index i, GNU_gama::Index j)
    {
      if (!this->is_solved) solve_me();
      return svd.q_bx(i, j) / this->sqrt_w(i);
    }
  
  void min_x()   {  svd.min_x(); }
  void min_x(GNU_gama::Index n, GNU_gama::Index x[]) { svd.min_x(n, x); }

  Float cond();
  
protected:
  
   void solve_me();
   
};

// ...................................................................

template <typename Float, typename Exc>
void OLSsvd<Float, Exc>::solve_me()
{
   using namespace GNU_gama; 
   using namespace std;

   if (this->is_solved) return;

   this->sqrt_w.reset(this->pb->dim());
   if (this->pw)
   { 
      for (GNU_gama::Index n = 1; n <= this->sqrt_w.dim(); n++)
         this->sqrt_w(n) = sqrt((*this->pw)(n));
      svd.reset(*this->pA, this->sqrt_w);                 // protected ==> public
      GNU_gama::Vec<Float, Exc> pbw = *this->pb;
      for (GNU_gama::Index i = 1; i <= this->pb->dim(); i++)
         pbw(i) *= this->sqrt_w(i);
      svd.solve(pbw, this->x);
   }
   else
   {
      svd.reset(*this->pA);
      for (GNU_gama::Index n = 1; n <= this->sqrt_w.dim(); n++)
         this->sqrt_w(n) = 1;
      svd.solve(*this->pb, this->x);
   }
   this->r  = *this->pA * this->x;
   this->r -= *this->pb;

   this->is_solved = true;
}

template <typename Float, typename Exc>
Float OLSsvd<Float, Exc>::cond()
{
  const GNU_gama::Vec<Float, Exc>& W = svd.SVD_W();

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









