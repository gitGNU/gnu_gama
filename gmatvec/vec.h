/*  
    C++ Matrix/Vector templates (gMatVec 0.9.14)
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
 *  $Id: vec.h,v 1.1 2001/12/07 11:59:46 cepek Exp $
 */

#ifndef gMatVec_Vec__h_
#define gMatVec_Vec__h_

#include <iostream>
#include <cstdarg>
#include <cmath>
#include <gmatvec/vecbase.h>

namespace gMatVec {
  

template <class Float=double, class Exc=Exception>
class Vec : public VecBase<Float, Exc> {

public:

  Vec() {}
  Vec(Index nsz) : VecBase<Float, Exc>(nsz) {}
  Vec(const VecBase<Float, Exc>& v) : VecBase<Float, Exc>(v) {}
  Vec(Index nsz, Float m11 ...) : VecBase<Float, Exc>(nsz)
    {
         using namespace std;
         MemRep<Float, Exc>::iterator p=begin();
         MemRep<Float, Exc>::iterator e=end();
         if (p == e)
            throw Exc(BadRank, "Vec::Vec(Index, Float ...)");
         *p = m11;  
         ++p;

         va_list  ap;
         va_start(ap, m11);
         while (p!=e)
            {
               *p = va_arg(ap, Float);
               ++p;
            }
         va_end(ap);
    }

  Vec operator*(Float f) const { 
    Vec t(dim()); mul(f, t); return t; 
  }
  Vec operator+(const Vec &x) const {
    Vec t(dim()); add(x, t); return t; 
  }
  Vec operator-(const Vec &x) const {
    Vec t(dim()); sub(x, t); return t; 
  }

  Vec& operator*=(Float f)      { mul(f, *this); return *this; }
  Vec& operator+=(const Vec &x) { add(x, *this); return *this; }
  Vec& operator-=(const Vec &x) { sub(x, *this); return *this; }

  MatVecBase<Float, Exc>::ListInitialiser operator=(Float x)
  {
    return list_init(x);
  }

};


template <class Float, class Exc>
inline Vec<Float, Exc> operator*(Float f, const Vec<Float, Exc>& V)
  {
    return V*f;
  }


template <class Float, class Exc>
Vec<Float, Exc> 
operator*(const MatBase<Float, Exc> &A, const Vec<Float, Exc> &b)
  {
    if (A.cols() != b.dim())
      throw Exc(BadRank, "Vec operator*(const MatBase&, const Vec&)");

    Vec<Float, Exc> t(A.rows());
    Float s;
    for (Index i=1; i<=A.rows(); i++)
      {
        s = 0;
        for (Index j=1; j<=A.cols(); j++) 
          s += A(i,j)*b(j);
        t(i) = s;
      }

    return t;
  }


template <class Float, class Exc>
Vec<Float, Exc> 
operator*(const Mat<Float, Exc> &A, const Vec<Float, Exc> &b)
  {
    if (A.cols() != b.dim())
      throw Exc(BadRank, "Vec operator*(const Mat&, const Vec&)");

    Vec<Float, Exc> t(A.rows());
    MemRep<Float, Exc>::iterator ti = t.begin();
    MemRep<Float, Exc>::const_iterator bb = b.begin();
    MemRep<Float, Exc>::const_iterator bi;
    MemRep<Float, Exc>::const_iterator ai = A.begin();
    Float s;
    for (Index i=1; i<=A.rows(); i++)
      {
        s = 0;
        bi = bb;
        for (Index j=1; j<=A.cols(); j++) 
          s += *ai++ * *bi++;
        *ti++ = s; 
      }

    return t;
  }


}   // namespace gMatVec

#endif
