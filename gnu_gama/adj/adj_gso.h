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
 *  $Id: adj_gso.h,v 1.1 2005/03/28 11:44:24 cepek Exp $
 */

#ifndef GNU_Gama_gnu_gama_gnugama_GaMa_OLS_gso_h
#define GNU_Gama_gnu_gama_gnugama_GaMa_OLS_gso_h

#include <gmatvec/gso.h>
#include <gnu_gama/adj/adj_base.h>
#include <cmath>

namespace GNU_gama {
  
template <typename Float, typename Exc>
class AdjGSO : public virtual AdjBase<Float, Exc> {

public:

  AdjGSO() {}
  AdjGSO(const GNU_gama::Mat<Float, Exc>& A, const GNU_gama::Vec<Float, Exc>& b)
    : AdjBase<Float, Exc>(A, b) {}
  AdjGSO(const GNU_gama::Mat<Float, Exc>& A, const GNU_gama::Vec<Float, Exc>& b,
         const GNU_gama::Vec<Float, Exc>& w) : AdjBase<Float, Exc>(A, b, w) {}
  
  void reset(const GNU_gama::Mat<Float, Exc>& A, 
             const GNU_gama::Vec<Float, Exc>& b)
    {
      AdjBase<Float, Exc>::reset(A, b);
    }
  void reset(const GNU_gama::Mat<Float, Exc>& A, 
             const GNU_gama::Vec<Float, Exc>& b,
             const GNU_gama::Vec<Float, Exc>& w)
    {
      AdjBase<Float, Exc>::reset(A, b, w);
    }
  
  const GNU_gama::Vec<Float, Exc>& solve(GNU_gama::Vec<Float, Exc>& x)
    {
      return x = AdjBase<Float, Exc>::solve();
    }
  const GNU_gama::Vec<Float, Exc>& solve() 
    { 
      return AdjBase<Float, Exc>::solve(); 
    }
  
  GNU_gama::Index defect() { return gso.defect(); }
  bool  lindep(GNU_gama::Index i) { return gso.lindep(i); }
  
  void  q_xx(GNU_gama::Mat<Float, Exc>& C) { AdjBase<Float, Exc>::q_xx(C); }
  Float q_xx(GNU_gama::Index i, GNU_gama::Index j);
  Float q_bb(GNU_gama::Index i, GNU_gama::Index j);
  Float q_bx(GNU_gama::Index i, GNU_gama::Index j);
  
  void min_x() { gso.min_x(); }
  void min_x(GNU_gama::Index n, GNU_gama::Index x[]) { gso.min_x(n, x); }
  
  
protected:
  
  void solve_me();
   
private:

   GNU_gama::Mat<Float, Exc> A_;
   GNU_gama::GSO<Float, Exc> gso;

   void init_gso_();
};

// ...................................................................


template <typename Float, typename Exc>
void AdjGSO<Float, Exc>::solve_me()
{
  if (this->is_solved) return;

  const GNU_gama::Index M = this->pA->rows();
  const GNU_gama::Index N = this->pA->cols();
  
  A_.reset(M+N, N+1);
  this->sqrt_w.reset(M);

  if (this->pw)
    {
      using namespace std;
      const GNU_gama::Vec<Float, Exc>& w_ = *this->pw;
      for (GNU_gama::Index i=1; i<=M; i++) this->sqrt_w(i) = sqrt(w_(i));
    }
  else
    {
      for (GNU_gama::Index i=1; i<=M; i++) this->sqrt_w(i) = 1;
    }

  const GNU_gama::Mat<Float, Exc>& A1 = *this->pA;
  const GNU_gama::Vec<Float, Exc>& b1 = *this->pb;

  {  // redundant curly braces needed by MS VC++ 
  for (GNU_gama::Index i=1; i<=M; i++)
    {
      A_(i, N+1) = -b1(i)*this->sqrt_w(i);
      for (GNU_gama::Index j=1; j<=N; j++) A_(i, j) = A1(i, j)*this->sqrt_w(i);
    }
  }
  {
  for (GNU_gama::Index i=1; i<=N; i++) 
    for (GNU_gama::Index j=1; j<=N+1; j++)
      A_(M+i, j) = (i==j) ? 1 : 0;
  }

  gso.reset(A_, M, N);
  gso.gso1(); 
  gso.gso2(); 

  this->x.reset(N);
  for (GNU_gama::Index i=1; i<=N; i++)
    this->x(i) = A_(M+i, N+1);

  this->r.reset(M);
  for (GNU_gama::Index j=1; j<=M; j++)
    this->r(j) = A_(j, N+1)/this->sqrt_w(j);
  
  this->is_solved = true; 
}


template <typename Float, typename Exc>
Float AdjGSO<Float, Exc>::q_xx(GNU_gama::Index i, GNU_gama::Index j)
  {
    if(!this->is_solved) solve_me();
    const GNU_gama::Index M = this->pA->rows();
    const GNU_gama::Index N = this->pA->cols();
    i += M;
    j += M;
    Float s = 0;                        
    for (GNU_gama::Index k=1; k<=N; k++) 
      s += A_(i,k)*A_(j,k);              // cov x_i x_j
    return s;
  }


template <typename Float, typename Exc>
Float AdjGSO<Float, Exc>::q_bb(GNU_gama::Index i, GNU_gama::Index j)
  {
    if(!this->is_solved) solve_me();

    const GNU_gama::Index N = this->pA->cols();
    Float s = 0;                        
    for (GNU_gama::Index k=1; k<=N; k++) 
      s += A_(i,k)*A_(j,k);              // cov b_i b_j
    return s/this->sqrt_w(i)/this->sqrt_w(j);
  }


template <typename Float, typename Exc>
Float AdjGSO<Float, Exc>::q_bx(GNU_gama::Index i, GNU_gama::Index j)
  {
    if(!this->is_solved) solve_me();

    const GNU_gama::Index M = this->pA->rows();
    const GNU_gama::Index N = this->pA->cols();
    j += M;
    Float s = 0;                        
    for (GNU_gama::Index k=1; k<=N; k++) 
      s += A_(i,k)*A_(j,k);              // cov b_i x_j
    return s/this->sqrt_w(i);
  }


}   // GaMaLib

#endif



















