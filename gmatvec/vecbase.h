/*  
    C++ Matrix/Vector templates (GNU Gama / gMatVec 0.9.23)
    Copyright (C) 1999  Ales Cepek <cepek@gnu.org>

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
 *  $Id: vecbase.h,v 1.13 2004/06/21 16:10:18 cepek Exp $
 *  http://www.gnu.org/software/gama/
 */

#ifndef gMatVec_VecBase__h_
#define gMatVec_VecBase__h_

#include <iostream>
#include <cstdarg>
#include <cmath>
#include <gmatvec/matvecbase.h>

namespace gMatVec {

  
template <typename Float=double, typename Exc=Exception>
class VecBase : public MatVecBase<Float, Exc> {

protected:

  VecBase() {}
  VecBase(Index nsz) :  MatVecBase<Float, Exc>(nsz) {}

public:

  typedef typename MatVecBase<Float, Exc>::iterator       iterator;
  typedef typename MatVecBase<Float, Exc>::const_iterator const_iterator;

  Index dim() const { return size(); }

  Float& operator()(Index n) { 
    Float* m = begin(); return m[--n]; 
  }
  Float  operator()(Index n) const { 
    const Float* m = begin(); return m[--n]; 
  } 

  void reset(Index n=0) { resize(n); }

  Float dot(const VecBase<Float, Exc> &B) const;

  Float norm_L1()   const;
  Float norm_L2()   const { return std::sqrt(dot(*this)); }
  Float norm_Linf() const;

};


template <typename Float, typename Exc>
Float VecBase<Float, Exc>::dot(const VecBase<Float, Exc> &B) const
  {
    if (dim() != B.dim())
      throw Exc(BadRank, "Float VecBase::dot(const VecBase&) const");
    
    const_iterator a = begin();
    const_iterator e = end();
    const_iterator b = B.begin();
    
    Float sum = 0;
    while (a != e) sum += *a++ * *b++;
    
    return sum;
  }


template <typename Float, typename Exc>
Float VecBase<Float, Exc>::norm_L1() const
  {
    const_iterator a = begin();
    const_iterator e = end();
    
    Float sum = 0;
    while (a != e) { sum += *a >= 0 ? *a : -(*a); ++a; }
    
    return sum;
  }


template <typename Float, typename Exc>
Float VecBase<Float, Exc>::norm_Linf() const
  {
    const_iterator a = begin();
    const_iterator e = end();
    
    Float norm = 0, x;
    while (a != e) 
      { 
        x = *a >= 0 ?  *a : -(*a); 
        ++a; 
        if (x > norm) norm = x;
      }
    
    return norm;
  }


template <typename Float, typename Exc>
std::istream& operator>>(std::istream& inp, VecBase<Float, Exc>& v)
  {
    Index size;
    inp >> size;
    
    if (size != v.dim())
      v.reset(size);
    
    typename MatVecBase<Float, Exc>::iterator b = v.begin();
    typename MatVecBase<Float, Exc>::iterator e = v.end();
    while (b != e)
      {
        inp >> *b;
        ++b;
      }
    
    return inp;
  }


}   // namespace gMatVec

#endif
