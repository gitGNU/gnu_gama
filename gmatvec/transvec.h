/*  
    C++ Matrix/Vector templates (GNU GaMa / gMatVec 0.9.17)
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
 *  $Id: transvec.h,v 1.4 2002/03/19 17:04:31 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec_TransVec__h_
#define gMatVec_TransVec__h_

#include <iostream>
#include <cstdarg>
#include <cmath>
#include <gmatvec/vec.h>

namespace gMatVec {
  

template <class Float=double, class Exc=Exception>
class TransVec : public VecBase<Float, Exc> {

public:

  typedef VecBase<Float, Exc>::iterator       iterator;
  typedef VecBase<Float, Exc>::const_iterator const_iterator;

  TransVec() {}
  TransVec(Index nsz) : VecBase<Float, Exc>(nsz) {}
  TransVec(const VecBase<Float, Exc>& v) : VecBase<Float, Exc>(v) {}
  
  TransVec operator*(Float f) const { 
      TransVec t(dim()); mul(f, t); return t; 
    }
  TransVec operator+(const TransVec &x) const {
    TransVec t(dim()); add(x, t); return t; 
  }
  TransVec operator-(const TransVec &x) const {
    TransVec t(dim()); sub(x, t); return t; 
  }

};


template <class Float, class Exc>
inline TransVec<Float, Exc> operator*(Float f, const TransVec<Float, Exc>& V)
  {
    return V*f;
  }


template <class Float, class Exc>
inline Float operator*(TransVec<Float, Exc> a, Vec<Float, Exc> b)
  {
    return a.dot(b);
  }


template <class Float, class Exc>
std::ostream& operator<<(std::ostream& out, const Vec<Float, Exc>& v)
  {
    const int fw = out.width();
    const int size = v.dim();
    out.width(fw);
    out << size << "\n\n";
    
    Vec<Float, Exc>::const_iterator b = v.begin();
    Vec<Float, Exc>::const_iterator e = v.end();
    while (b != e)
      {
        out.width(fw);
        out << *b;
        ++b;
        out << '\n';
      }
    
    out << '\n';
    return out;
  }


template <class Float, class Exc>
std::ostream& operator<<(std::ostream& out, const TransVec<Float, Exc>& v)
  {
    const int fw = out.width();
    const int size = v.dim();
    out.width(fw);
    out << size << "  ";
    
    TransVec<Float, Exc>::const_iterator b = v.begin();
    TransVec<Float, Exc>::const_iterator e = v.end();
    while (b != e)
      {
        out.width(fw);
        out << *b;
        ++b;
        out << ' ';
      }
    
    out << '\n';
    return out;
  }


template <class Float, class Exc>
inline TransVec<Float, Exc> trans(const Vec<Float, Exc>& v) 
{ 
  TransVec<Float, Exc> T(v); return T; 
}


template <class Float, class Exc>
inline Vec<Float, Exc> trans(const TransVec<Float, Exc>& v) 
{ 
  Vec<Float, Exc> T(v); return T; 
}


template <class Float, class Exc>
TransVec<Float, Exc> 
operator*(const TransVec<Float, Exc> &b, const MatBase<Float, Exc> &A)
  {
    if (b.dim() != A.rows())
      throw Exc(BadRank, "TransVec operator*(const TransVec&, const MatBase&)");

    TransVec<Float, Exc> t(A.cols());
    Float s;
    for (Index j=1; j<=A.cols(); j++)
      {
        s = 0;
        for (Index i=1; i<=A.cols(); i++) 
          s += b(i)*A(i,j);
        t(j) = s;
      }

    return t;
  }


template <class Float, class Exc>
TransVec<Float, Exc> 
operator*(const TransVec<Float, Exc> &b, const Mat<Float, Exc> &A)
  {
    if (b.dim() != A.rows())
      throw Exc(BadRank, "TransVec operator*(const TransVec&, const Mat&)");

    TransVec<Float, Exc> t(A.cols());
    TransVec<Float, Exc>::iterator ti =t.begin();
    TransVec<Float, Exc>::const_iterator bb = b.begin();
    TransVec<Float, Exc>::const_iterator bi;
    Mat<Float, Exc>::const_iterator aj = A.begin();
    Mat<Float, Exc>::const_iterator ai;
    const Index a_cols = A.cols();
    const Index a_rows = A.rows();
    Float s;
    for (Index j=1; j<=a_cols; j++)
      {
        s = 0;
        bi = bb;
        ai = aj;
        for (Index i=0; i<a_rows; i++, ai += a_cols) 
          s += *bi++ * *ai;
        *ti++ = s;
        aj++;
      }

    return t;
  }


}   // namespace gMatVec

#endif








